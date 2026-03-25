from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from pyNJOY2GENDF.gendf_parser import GENDFparser
from pyNJOY2GENDF.interpolation import interp1, interpolate_sig0, interpolate_temperature, IsoSig0Interpolated


@dataclass
class MaterialSpec:
    """Specification for a single material to process."""
    name: str
    nuclides: list[str]  # GENDF filenames (e.g., ["U235.gendf", "U238.gendf"])
    number_densities: list[float]  # atoms/barn-cm
    temperature: float  # Kelvin
    gendf_directory: Path


class MatXS:
    """Extract macroscopic cross-sections from GENDF files for materials.
    
    Primary usage: instantiate with MaterialSpec objects and call process().
    For legacy MONTE input files, use the from_monte_input() classmethod.
    """
    
    def __init__(self, materials: list[MaterialSpec], printmicro: bool = False):
        """Initialize with material specifications.
        
        Args:
            materials: List of MaterialSpec objects defining materials to process
            printmicro: If True, write microscopic XS per isotope alongside macro XS
        """
        self.materials = materials
        self.printmicro = printmicro
        self._iso_cache: dict[tuple[Path, str], dict] = {}
    
    @classmethod
    def from_monte_input(cls, input_path: str | Path, printmicro: bool = False) -> MatXS:
        """Create from a MONTE input deck file.
        
        Args:
            input_path: Path to MONTE input file with XSP/SIGNAL/MAT cards
            printmicro: If True, write microscopic XS per isotope
            
        Returns:
            MatXS instance ready to process()
        """
        materials = cls._parse_monte_input(Path(input_path))
        return cls(materials, printmicro)
    
    @staticmethod
    def _parse_monte_input(path: Path) -> list[MaterialSpec]:
        """Parse MONTE input file format."""
        gendf_path: Path | None = None
        signals: dict[str, float] = {}  # signal_name -> value
        mat_components: dict[str, list[tuple[str, float, str]]] = {}  # matname -> [(iso, numden, temp_signal)]

        with path.open("r") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                first = line[0]
                if first == "$":
                    break
                if first == "*":
                    continue

                parts = line.split()
                kw = parts[0].upper()

                if kw == "XSP":
                    p = parts[1].strip("\"'")
                    gendf_path = Path(p)

                elif kw in ("SIG", "SIGNAL"):
                    name = parts[1]
                    kind = parts[2].upper()
                    if kind == "CONSTANT":
                        value = float(parts[3])
                    else:
                        raise ValueError(f"Unsupported signal type: {kind}")
                    signals[name] = value

                elif kw == "MAT":
                    matname = parts[1]
                    isoname = parts[2]
                    numden = float(parts[3])
                    tsig = parts[4]

                    if matname not in mat_components:
                        mat_components[matname] = []
                    mat_components[matname].append((isoname, numden, tsig))

        if gendf_path is None:
            raise ValueError("Input deck missing XSP path")

        # Convert to MaterialSpec objects
        material_specs: list[MaterialSpec] = []
        for matname, components in mat_components.items():
            nuclides = [c[0] for c in components]
            numdensities = [c[1] for c in components]
            temp_signal_name = components[0][2]  # assume all use same temp
            temperature = signals.get(temp_signal_name)
            if temperature is None:
                raise ValueError(f"Temperature signal {temp_signal_name} not defined")
            
            material_specs.append(
                MaterialSpec(
                    name=matname,
                    nuclides=nuclides,
                    number_densities=numdensities,
                    temperature=temperature,
                    gendf_directory=gendf_path,
                )
            )
        
        return material_specs
    
    def _load_isotope(self, gendf_directory: Path, filename: str) -> dict:
        """Load all data from a GENDF file (cached)."""
        # Check cache
        cache_key = (gendf_directory, filename)
        if cache_key in self._iso_cache:
            return self._iso_cache[cache_key]
        
        p = gendf_directory / filename
        parser = GENDFparser(p)

        ng = parser.get_ng()
        nsig0 = parser.get_nsig0()
        sig0_grid = parser.get_sig0_grid()
        temps = parser.get_temperatures()

        # MF=5 MT=18 chi (may be absent)
        chi = parser.extract_mf5_mt18_chi()

        # MF=3 MT=452 nubar (may be absent)
        nubar, _ = parser.extract_mf3_mt(mt=452, itemp=0, nl=0)
        nubar_vec = nubar[:, 0] if nubar.shape[1] else np.zeros((ng,), dtype=float)

        # Kerma: MF=3 MT=301, take xs part for first temperature (itemp=0)
        kerma_xs, _ = parser.extract_mf3_mt(mt=301, itemp=0, nl=0)

        # Temperature-dependent sigt/sigf/flux from MF=3 MT=1 and MT=18
        ntemp = max(1, len(temps))
        sigt = np.zeros((ntemp, ng, nsig0), dtype=float)
        sigf = np.zeros((ntemp, ng, nsig0), dtype=float)
        flux = np.zeros((ntemp, ng, nsig0), dtype=float)

        for it in range(ntemp):
            xs_t, fl_t = parser.extract_mf3_mt(mt=1, itemp=it, nl=0)
            sigt[it] = xs_t
            flux[it] = fl_t

            xs_f, _ = parser.extract_mf3_mt(mt=18, itemp=it, nl=0)
            sigf[it] = xs_f

        # Elastic MF=6 MT=2 for all Legendre orders and temperatures
        # Note: In some GENDF files, elastic scatter may only be present at first temperature
        # In that case, reuse the pattern for all temperatures
        f_e = np.zeros((0,), dtype=int)
        t_e = np.zeros((0,), dtype=int)
        sige = None

        # Determine maximum nl by trying to extract increasing orders
        nl_el = 0
        for test_nl in range(20):
            try:
                f_test, t_test, sig_test = parser.extract_mf6_scatter(mt=2, itemp=0, nl=test_nl)
                if len(f_test) > 0:
                    nl_el = test_nl + 1
                else:
                    break
            except Exception:
                break
        
        if nl_el == 0:
            # No elastic scatter data
            f_e = np.zeros((0,), dtype=int)
            t_e = np.zeros((0,), dtype=int)
            sige = np.zeros((ntemp, 1, 0, nsig0), dtype=float)
            nl_el = 1
        else:
            try:
                # Extract first Legendre order to get sparsity pattern
                f_tmp, t_tmp, sig_tmp = parser.extract_mf6_scatter(mt=2, itemp=0, nl=0)
                nonze = len(f_tmp)
                f_e = f_tmp
                t_e = t_tmp
                sige = np.zeros((ntemp, nl_el, nonze, nsig0), dtype=float)
                
                # Extract all Legendre orders for all temperatures
                for ell in range(nl_el):
                    for it in range(ntemp):
                        _, _, sig_it = parser.extract_mf6_scatter(mt=2, itemp=it, nl=ell)
                        if len(sig_it) == 0:
                            # No elastic scatter at this temperature, reuse first temperature
                            _, _, sig_it = parser.extract_mf6_scatter(mt=2, itemp=0, nl=ell)
                        if sig_it.shape[0] != nonze:
                            raise ValueError(f"Elastic nonzero pattern changed with Legendre order {ell+1} in {filename}")
                        sige[it, ell] = sig_it
            except Exception:
                # If elastic scatter extraction fails, proceed with empty arrays
                nonze = 0
                f_e = np.zeros((0,), dtype=int)
                t_e = np.zeros((0,), dtype=int)
                sige = np.zeros((ntemp, 1, 0, nsig0), dtype=float)
                nl_el = 1

        # Inelastic MF=6 MT=51..91 (temperature-independent, extract at itemp=0 only)
        # Determine max nl for inelastic (should match elastic, but check first MT)
        nl_inel = nl_el  # Default to same as elastic
        for mt in range(51, 92):
            try:
                f_test, _, _ = parser.extract_mf6_scatter(mt=mt, itemp=0, nl=0, strict=False)
                if len(f_test) > 0:
                    # Found inelastic data, determine its max nl
                    for test_nl in range(20):
                        try:
                            f_t, _, _ = parser.extract_mf6_scatter(mt=mt, itemp=0, nl=test_nl, strict=False)
                            if len(f_t) == 0:
                                nl_inel = test_nl
                                break
                        except:
                            nl_inel = test_nl
                            break
                    break
            except:
                continue
        
        # Extract all inelastic MTs for all Legendre orders
        f_i_list: list[int] = []
        t_i_list: list[int] = []
        sigi_per_nl: list[list[np.ndarray]] = [[] for _ in range(max(nl_inel, 1))]
        
        for mt in range(51, 92):
            for ell in range(max(nl_inel, 1)):
                try:
                    f_mt, t_mt, sig_mt = parser.extract_mf6_scatter(mt=mt, itemp=0, nl=ell, strict=False)
                    if len(f_mt) == 0:
                        continue
                    if ell == 0:  # Only collect f,t once
                        f_i_list.extend(f_mt)
                        t_i_list.extend(t_mt)
                    sigi_per_nl[ell].append(sig_mt)
                except:
                    continue

        f_i = np.asarray(f_i_list, dtype=int)
        t_i = np.asarray(t_i_list, dtype=int)
        
        # Stack into (nl, nonzero, nsig0) array
        if any(sigi_per_nl):
            sigi = np.zeros((max(nl_inel, 1), len(f_i_list), nsig0), dtype=float)
            for ell in range(max(nl_inel, 1)):
                if sigi_per_nl[ell]:
                    sigi[ell] = np.vstack(sigi_per_nl[ell])
        else:
            sigi = np.zeros((1, 0, nsig0), dtype=float)

        result = {
            "isoname": filename,
            "temps": temps,
            "sig0_grid": sig0_grid,
            "sigt": sigt,
            "sigf": sigf,
            "kerma": kerma_xs,
            "flux": flux,
            "nubar": nubar_vec,
            "chi": chi,
            "f_e": f_e,
            "t_e": t_e,
            "sige": sige,
            "f_i": f_i,
            "t_i": t_i,
            "sigi": sigi,
        }
        
        self._iso_cache[cache_key] = result
        return result
    
    @staticmethod
    def _sigma0_iteration(
        iso_t_list, numdens: np.ndarray, tol: float = 1.0e-4, maxiter: int = 1000
    ) -> tuple[np.ndarray, list[IsoSig0Interpolated]]:
        """Perform sigma0 iteration for multi-isotope materials."""
        # iso_t_list: list of IsoTempInterpolated, each has sigt (ng,nsig0)
        niso = len(iso_t_list)
        ng = iso_t_list[0].ng

        sig0 = np.ones((niso, ng), dtype=float)
        sig0_out = np.ones((niso, ng), dtype=float)

        iso_s_list: list[IsoSig0Interpolated] = []

        if niso == 1:
            iso_s_list = [interpolate_sig0(iso_t_list[0], sig0[0])]
            return sig0, iso_s_list

        err = 1.0
        it = 0
        
        print(f"=== DEBUG MatXS: Starting sig0 iteration ===")
        print(f"    niso={niso}, ng={ng}")
        print(f"    numdens={numdens}")
        
        while err > tol:
            iso_s_list = [interpolate_sig0(iso_t_list[i], sig0[i]) for i in range(niso)]

            # DEBUG: Print initial interpolated total XS
            if it == 0:
                print(f"    Initial interpolated total XS:")
                for i in range(niso):
                    print(f"      Iso {i+1}: {iso_s_list[i].sigt}")
            
            err = 0.0
            for g in range(ng):
                for i in range(niso):
                    s = 0.0
                    for j in range(niso):
                        if j == i:
                            continue
                        s += numdens[j] * iso_s_list[j].sigt[g]
                    s = s / numdens[i]
                    s = min(s, 1.0e10)
                    s = max(s, 1.0)
                    sig0_out[i, g] = s

                err = max(err, float(np.max(np.abs(sig0[:, g] - sig0_out[:, g]))))

            sig0[:, :] = sig0_out[:, :]
            it += 1
            
            # DEBUG: Print convergence
            if it <= 3 or err <= tol:
                print(f"    Iteration {it}, error={err:.6e}")
                if it <= 3:
                    for i in range(niso):
                        print(f"      Iso {i+1} sig0: {sig0[i, :]}")
            
            if it > maxiter:
                raise RuntimeError("Too many sigma0 iterations")

        print(f"=== DEBUG MatXS: Final sig0 after convergence ===")
        for i in range(niso):
            print(f"    Iso {i+1}:")
            for g in range(ng):
                print(f"      Group {g+1}: {sig0[i,g]:.8e}")
        print()
        
        iso_s_list = [interpolate_sig0(iso_t_list[i], sig0[i]) for i in range(niso)]
        return sig0, iso_s_list
    
    def _compute_material(self, mat: MaterialSpec) -> dict:
        """Compute macroscopic cross-sections for one material.
        
        Returns:
            Dictionary containing material data with structure:
            {
                'name': str,
                'temperature': float,
                'nuclides': list[str],
                'number_densities': ndarray,
                'sig0': ndarray (niso, ng),
                'macro': {
                    'total': ndarray (ng,),
                    'fission': ndarray (ng,),
                    'nu_fission': ndarray (ng,),
                    'chi': ndarray (ng,),
                    'kerma': ndarray (ng,),
                    'absorption': ndarray (ng,),
                    'scatter': ndarray (nl, ng, ng)
                },
                'micro': {  # Only if printmicro=True
                    '<nuclide_name>': {
                        'total': ndarray,
                        'fission': ndarray,
                        'kerma': ndarray,
                        'nu': ndarray,
                        'chi': ndarray,
                        'flux': ndarray (optional)
                    }, ...
                }
            }
        """
        # Build iso_t for each component
        iso_names = mat.nuclides
        numdens = np.array(mat.number_densities, dtype=float)

        iso_t_list = []
        for iso_filename in mat.nuclides:
            iso_raw = self._load_isotope(mat.gendf_directory, iso_filename)

            iso_t = interpolate_temperature(
                isoname=iso_filename,
                temps=iso_raw["temps"],
                sig0_grid=iso_raw["sig0_grid"],
                sigt=iso_raw["sigt"],
                sigf=iso_raw["sigf"],
                kerma=iso_raw["kerma"],
                flux=iso_raw["flux"],
                nubar=iso_raw["nubar"],
                chi=iso_raw["chi"],
                temp_query=mat.temperature,
                f_e=iso_raw["f_e"],
                t_e=iso_raw["t_e"],
                sige=iso_raw["sige"],
                f_i=iso_raw["f_i"],
                t_i=iso_raw["t_i"],
                sigi=iso_raw["sigi"],
            )
            iso_t_list.append(iso_t)

        sig0, iso_s_list = self._sigma0_iteration(iso_t_list, numdens)

        ng = iso_s_list[0].ng
        nl_max = max(int(iso.nl) for iso in iso_s_list)

        SigT = np.zeros((ng,), dtype=float)
        SigF = np.zeros((ng,), dtype=float)
        SigP = np.zeros((ng,), dtype=float)
        Chi = np.zeros((ng,), dtype=float)
        Kerma = np.zeros((ng,), dtype=float)

        for i, iso in enumerate(iso_s_list):
            SigT += numdens[i] * iso.sigt
            SigF += numdens[i] * iso.sigf
            SigP += numdens[i] * iso.nubar * iso.sigf
            Kerma += numdens[i] * iso.kerma

        # Chi (legacy behavior): weighted by number density where sigf>0
        for g in range(ng):
            s = 0.0
            for i, iso in enumerate(iso_s_list):
                if iso.sigf[g] > 0:
                    s += numdens[i] * iso.chi[g]
            Chi[g] = s
        if Chi.sum() > 0:
            Chi /= Chi.sum()

        SigS = np.zeros((nl_max, ng, ng), dtype=float)

        for i, iso in enumerate(iso_s_list):
            for ell in range(min(nl_max, iso.nl)):
                # elastic
                for j in range(len(iso.f_e)):
                    f = int(iso.f_e[j]) - 1
                    t = int(iso.t_e[j]) - 1
                    SigS[ell, f, t] += numdens[i] * iso.sige[ell, j]
                # inelastic
                for j in range(len(iso.f_i)):
                    f = int(iso.f_i[j]) - 1
                    t = int(iso.t_i[j]) - 1
                    SigS[ell, f, t] += numdens[i] * iso.sigi[ell, j]

        SigST = SigS[0].sum(axis=1)
        SigA = SigT - SigST

        # Build result dictionary
        result = {
            'name': mat.name,
            'temperature': mat.temperature,
            'nuclides': list(mat.nuclides),
            'number_densities': numdens,
            'sig0': sig0,
            'macro': {
                'total': SigT,
                'fission': SigF,
                'nu_fission': SigP,
                'chi': Chi,
                'kerma': Kerma,
                'absorption': SigA,
                'scatter': SigS
            }
        }

        if self.printmicro:
            micro = {}
            for i, iso in enumerate(iso_s_list):
                iso_data = {
                    'total': iso.sigt,
                    'fission': iso.sigf,
                    'kerma': iso.kerma,
                    'nu': iso.nubar,
                    'chi': iso.chi
                }
                if iso.flux is not None:
                    iso_data['flux'] = iso.flux
                micro[iso.isoname] = iso_data
            result['micro'] = micro

        return result
    
    def process(self) -> dict[str, dict]:
        """Process all materials and return results as dictionaries.
        
        Returns:
            Dictionary mapping material names to their computed data:
            {
                '<material_name>': {
                    'name': str,
                    'temperature': float,
                    'nuclides': list[str],
                    'number_densities': ndarray,
                    'sig0': ndarray,
                    'macro': {...},
                    'micro': {...}  # Only if printmicro=True
                }, ...
            }
        """
        results = {}
        
        # Process each material
        for mat in self.materials:
            # Load isotopes (cached)
            for iso_filename in mat.nuclides:
                print(f"Processing GENDF file: {mat.gendf_directory / iso_filename}")
                self._load_isotope(mat.gendf_directory, iso_filename)
            
            print(f"Processing material sig0 and temp interpolation: {mat.name}")
            results[mat.name] = self._compute_material(mat)
        
        return results
    
    def process_to_csv(self, output_directory: Optional[str | Path] = None) -> dict[str, dict]:
        """Process all materials and write output CSV files.
        
        Args:
            output_directory: Where to write material_xs/ outputs. 
                            If None, writes to current working directory.
        
        Returns:
            Dictionary mapping material names to their computed data (same as process())
        """
        # Get results
        results = self.process()
        
        # Setup output directory
        if output_directory is None:
            outdir = Path.cwd()
        else:
            outdir = Path(output_directory)
        
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "material_xs").mkdir(exist_ok=True)
        
        # Write material list
        with (outdir / "material_xs" / "material_list").open("w") as f:
            f.write("Material\n")
            for mat_name in results.keys():
                f.write(f"{mat_name}\n")
        
        # Write CSV files for each material
        for mat_name, mat_data in results.items():
            base = outdir / "material_xs"
            iso_names = mat_data['nuclides']
            _ensure_material_dirs(base, mat_name, iso_names, self.printmicro)
            
            # Write macro files
            macro = base / mat_name / "macro"
            _write_csv_1d(macro / "total.csv", mat_data['macro']['total'])
            _write_csv_1d(macro / "fission.csv", mat_data['macro']['fission'])
            _write_csv_1d(macro / "nu_fission.csv", mat_data['macro']['nu_fission'])
            _write_csv_1d(macro / "chi.csv", mat_data['macro']['chi'])
            _write_csv_1d(macro / "kerma.csv", mat_data['macro']['kerma'])
            _write_csv_1d(macro / "absorption.csv", mat_data['macro']['absorption'])
            _write_sparse_scatter(macro / "scatter.csv", mat_data['macro']['scatter'])
            
            # Write isolist
            with (base / mat_name / "isolist.csv").open("w") as f:
                f.write("Iso,Temp,Conc\n")
                for iso_filename, nd in zip(mat_data['nuclides'], mat_data['number_densities']):
                    f.write(f"{iso_filename},{mat_data['temperature']:.14e},{nd:.14e}\n")
            
            # Write sig0 values per isotope
            sig0 = mat_data['sig0']
            ng = sig0.shape[1]
            with (base / mat_name / "sig0.csv").open("w") as f:
                f.write("Iso,GroupIndex,Value\n")
                for i in range(len(mat_data['nuclides'])):
                    for g in range(ng):
                        f.write(f"{i+1},{g+1},{sig0[i,g]:.14e}\n")
            
            # Write microscopic files if requested
            if self.printmicro and 'micro' in mat_data:
                for iso_name, iso_data in mat_data['micro'].items():
                    iso_dir = base / mat_name / iso_name
                    _write_csv_1d(iso_dir / "total.csv", iso_data['total'])
                    _write_csv_1d(iso_dir / "fission.csv", iso_data['fission'])
                    _write_csv_1d(iso_dir / "kerma.csv", iso_data['kerma'])
                    _write_csv_1d(iso_dir / "nu.csv", iso_data['nu'])
                    _write_csv_1d(iso_dir / "chi.csv", iso_data['chi'])
                    if 'flux' in iso_data:
                        _write_csv_1d(iso_dir / "scalarflux.csv", iso_data['flux'])
        
        return results


# Helper functions
def _write_csv_1d(path: Path, arr: np.ndarray) -> None:
    """Write 1D array to CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("GroupIndex,Value\n")
        for i, v in enumerate(arr, start=1):
            f.write(f"{i},{v:.14e}\n")


def _write_sparse_scatter(path: Path, SigS: np.ndarray, threshold: float = 1.0e-15) -> None:
    """Write sparse scattering matrix to CSV."""
    # SigS: (nl, ng, ng)
    nl, ng, _ = SigS.shape
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("LegendreIndex,FromIndex,ToIndex,Value\n")
        for ell in range(nl):
            A = SigS[ell]
            # sparse scan
            for gfrom in range(ng):
                row = A[gfrom]
                nz = np.where(np.abs(row) > threshold)[0]
                for gto in nz:
                    f.write(f"{ell+1},{gfrom+1},{gto+1},{row[gto]:.14e}\n")


def _ensure_material_dirs(base: Path, matname: str, iso_names: list[str], printmicro: bool) -> None:
    """Create directory structure for material outputs."""
    (base / matname / "macro").mkdir(parents=True, exist_ok=True)
    if printmicro:
        for iso in iso_names:
            (base / matname / iso).mkdir(parents=True, exist_ok=True)


def main() -> None:
    """CLI entry point."""
    ap = argparse.ArgumentParser(description="Extract macroscopic cross-sections from GENDF files")
    ap.add_argument("--input", default="input", help="Path to MONTE input deck (default: ./input)")
    ap.add_argument("--printmicro", default="false", choices=["true", "false"], help="Write microscopic XS")
    args = ap.parse_args()

    printmicro = args.printmicro.lower() == "true"
    inp = Path(args.input)
    
    # Create from MONTE input file
    processor = MatXS.from_monte_input(inp, printmicro=printmicro)
    
    # Process and write to directory containing input file
    processor.process_to_csv(output_directory=inp.resolve().parent)


if __name__ == "__main__":
    main()
