# pyNJOY2GENDF
A Python package to automate NJOY execution for generating multigroup cross-section libraries (GENDF format) and extracting macroscopic cross-sections for reactor physics applications.

## Features

- **NJOY Input Generation**: Automated creation of NJOY input files for RECONR, BROADR, UNRESR/PURR, HEATR, and GROUPR modules
- **GENDF Parsing**: Robust parser for reading multigroup cross-section data from GENDF files
- **Macroscopic XS Extraction**: Calculate material-level cross-sections with temperature and sigma0 interpolation
- **Parallel Execution**: Multi-threaded NJOY processing for large isotope libraries
- **Flexible API**: Direct material specification or legacy MONTE input file support

## Installation 

Install this python package in edit mode with the command in this folder:

```bash
pip3 install -e .
```

## Quick Start

### Generating GENDF Files

```python
import pyNJOY2GENDF.njoy2gendf

# Create NJOY input generator
wr = pyNJOY2GENDF.njoy2gendf.njoy2gendf(
    savedir='gendf_library',
    path2endf='path/to/endf/files',
    njoy_exec='njoy21',
    nuclide_list=['U235', 'U238', 'Pu239'],
    filename_convention='ENDF-B'
)

# Configure and write NJOY inputs
wr.write_njoy_gendf_inputs(
    energybin_edges=[1e-5, 1, 100, 1000, 2e6],
    library_name='ENDF-B-VIII.0',
    temperatures=[300, 600, 900],
    sig0s=[1e10, 1e5, 1e4, 1000, 100, 10, 1],
    groupr_iwt=8  # Thermal reactor spectrum
)

# Run NJOY (parallel execution)
wr.run_njoy_inputs(nthreads=4)
```

### Extracting Macroscopic Cross-Sections

```python
from pyNJOY2GENDF.MatXS import MatXS, MaterialSpec
from pathlib import Path

# Define material composition
fuel = MaterialSpec(
    name='UO2_FUEL',
    nuclides=['U235.gendf', 'U238.gendf', 'O16.gendf'],
    number_densities=[4.5e-3, 2.0e-2, 4.2e-2],  # atoms/barn-cm
    temperature=900.0,  # Kelvin
    gendf_directory=Path('gendf_library/GENDF')
)

# Create processor and extract macroscopic XS
processor = MatXS([fuel], printmicro=False)
processor.process(output_directory='macro_xs')
```

## API Reference

### MatXS Class

Extract macroscopic cross-sections from GENDF files with temperature and sigma0 interpolation.

#### Public Methods

**`__init__(materials: List[MaterialSpec], printmicro: bool = False)`**
- Initialize with material specifications
- `materials`: List of MaterialSpec objects defining materials to process
- `printmicro`: If True, write microscopic XS per isotope alongside macro XS

**`from_monte_input(input_path: str | Path, printmicro: bool = False) -> MatXS`** (classmethod)
- Create from a MONTE input deck file
- `input_path`: Path to MONTE input file with XSP/SIGNAL/MAT cards
- `printmicro`: If True, write microscopic XS per isotope
- Returns: MatXS instance ready to process()

**`process(output_directory: Optional[str | Path] = None) -> None`**
- Process all materials and write output files
- `output_directory`: Where to write material_xs/ outputs (default: current directory)

#### Private Methods

**`_parse_monte_input(path: Path) -> List[MaterialSpec]`** (static)
- Parse MONTE input file format
- Returns list of MaterialSpec objects

**`_load_isotope(gendf_directory: Path, filename: str) -> dict`**
- Load and cache all data from a GENDF file
- Returns dictionary with cross-section arrays and metadata

**`_sigma0_iteration(iso_t_list, numdens: np.ndarray, tol: float, maxiter: int)`** (static)
- Perform sigma0 iteration for multi-isotope materials
- Returns: (sig0 array, list of interpolated isotope data)

**`_compute_material(mat: MaterialSpec, outdir: Path) -> None`**
- Compute and write macroscopic cross-sections for one material
- Generates CSV files for total, fission, scattering, etc.

#### Helper Functions (Module-level)

**`_write_csv_1d(path: Path, arr: np.ndarray) -> None`**
- Write 1D cross-section array to CSV file
- Format: "GroupIndex,Value"

**`_write_sparse_scatter(path: Path, SigS: np.ndarray, threshold: float) -> None`**
- Write sparse scattering matrix to CSV
- Format: "LegendreIndex,FromIndex,ToIndex,Value"

**`_ensure_material_dirs(base: Path, matname: str, iso_names: List[str], printmicro: bool) -> None`**
- Create directory structure for material outputs

### MaterialSpec Dataclass

Specification for a single material to process.

**Attributes:**
- `name: str` - Material name
- `nuclides: List[str]` - GENDF filenames (e.g., ["U235.gendf", "U238.gendf"])
- `number_densities: List[float]` - Atom densities in atoms/barn-cm
- `temperature: float` - Material temperature in Kelvin
- `gendf_directory: Path` - Path to directory containing GENDF files

## Output Files

The `MatXS.process()` method generates the following directory structure:

```
output_directory/
└── material_xs/
    ├── material_list           # List of all materials
    └── <material_name>/
        ├── isolist.csv         # Isotope list with temperatures and densities
        ├── sig0.csv            # Background cross-sections per isotope and group
        └── macro/
            ├── total.csv       # Total cross-section
            ├── fission.csv     # Fission cross-section
            ├── nu_fission.csv  # Nu*fission cross-section
            ├── chi.csv         # Fission spectrum
            ├── absorption.csv  # Absorption cross-section
            ├── kerma.csv       # Kinetic energy release per mass
            └── scatter.csv     # Scattering matrix (sparse format)
```

## Examples

See the `examples/` directory for complete working examples:
- `endfb8_example/` - ENDF/B-VIII.0 library generation
- `jeff4_example/` - JEFF-4.0 library generation
- `matxs_example/` - Macroscopic cross-section extraction

## Dependencies

- numpy >= 2.3.3
- tabulate
- Python >= 3.13.9
- NJOY21 executable (for GENDF generation)

## Citation

If you use this package in your research, please cite the relevant nuclear data libraries and NJOY code.

