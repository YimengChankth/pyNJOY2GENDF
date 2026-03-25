from __future__ import annotations

from dataclasses import dataclass
import math
import re
from pathlib import Path
from typing import Iterable, Iterator, Optional

import numpy as np


_ENDF_EXP_RE = re.compile(r"^([+-]?\d*\.?\d*)([+-]\d+)$")


def _parse_endf_number(field: str) -> Optional[float]:
    """Parse an 11-character ENDF numeric field.

    Handles ENDF's implicit exponent format, e.g. ' 1.234567-3' -> 1.234567E-3.
    Returns None for blank fields.
    """

    s = field.strip()
    if not s:
        return None

    # Already standard float
    if "e" in s.lower():
        return float(s)

    # Integers (ENDF also uses integer fields in the same 11-char slots)
    if "." not in s and ("+" not in s[1:] and "-" not in s[1:]):
        return float(int(s))

    # Implicit exponent format: split mantissa and trailing exponent sign+digits
    m = _ENDF_EXP_RE.match(s)
    if m and m.group(2):
        mant, exp = m.group(1), m.group(2)
        if mant in ("", "+", "-"):
            # Defensive: malformed
            return float(s)
        return float(f"{mant}E{exp}")

    return float(s)


@dataclass(frozen=True)
class ENDFRecord:
    w: tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]
    mf: int
    mt: int
    ns: int  # line number within section, ENDF columns 76-80

    @property
    def c1(self) -> float:
        return float(self.w[0] or 0.0)

    @property
    def c2(self) -> float:
        return float(self.w[1] or 0.0)

    @property
    def l1(self) -> int:
        return int(self.w[2] or 0)

    @property
    def l2(self) -> int:
        return int(self.w[3] or 0)

    @property
    def n1(self) -> int:
        return int(self.w[4] or 0)

    @property
    def n2(self) -> int:
        return int(self.w[5] or 0)


class GENDFparser:
    """Lightweight parser for NJOY GENDF (ENDF-6 formatted) groupwise files.

    This is intentionally pragmatic (not a full ENDF-6 implementation): it focuses on
    the record patterns used by NJOY's GROUPR GENDF outputs.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.title: str = ""
        self.records: list[ENDFRecord] = []

        self._load()

    def _load(self) -> None:
        with self.path.open("r") as f:
            first = f.readline()
            if not first:
                raise ValueError(f"Empty GENDF file: {self.path}")
            self.title = first.rstrip("\n")

            for raw in f:
                line = raw.rstrip("\n")
                if not line:
                    continue
                if len(line) < 80:
                    line = line.ljust(80)

                # End-of-file sentinel used in your Fortran: line(69:70) == '-1'
                # (1-based). That's [68:70) 0-based.
                if line[68:70].strip() == "-1":
                    break

                fields = [line[i * 11 : (i + 1) * 11] for i in range(6)]
                w = tuple(_parse_endf_number(x) for x in fields)  # type: ignore[assignment]

                mf = int(line[70:72])
                mt = int(line[72:75])
                ns = int(line[75:80])
                self.records.append(ENDFRecord(w=w, mf=mf, mt=mt, ns=ns))

        if not self.records:
            raise ValueError(f"No records parsed from {self.path}")

    def _section_starts(self, mf: int, mt: int) -> list[int]:
        return [
            i
            for i, r in enumerate(self.records)
            if r.mf == mf and r.mt == mt and r.ns == 1
        ]

    def get_ng(self) -> int:
        # Matches the legacy Fortran assumption: cards(2,3) (1-based).
        if len(self.records) < 2:
            raise ValueError("GENDF too short to infer ng")
        return int(self.records[1].l1)

    def get_nsig0(self) -> int:
        # Matches Fortran: cards(1,4) (1-based): first record after title, 4th word.
        return int(self.records[0].l2)

    def get_sig0_grid(self) -> np.ndarray:
        nsig0 = self.get_nsig0()
        # Matches Fortran: read nsig0+1 words starting at cards row 3 (1-based) => records[2]
        start = 2
        words = self.read_words(start, nsig0 + 1)
        return words[1:].copy()

    def get_temperatures(self) -> np.ndarray:
        starts = self._section_starts(mf=1, mt=451)
        temps = []
        for s in starts:
            if s + 1 >= len(self.records):
                continue
            temps.append(self.records[s + 1].c1)
        return np.array(temps, dtype=float)

    def read_words(self, start_record_index: int, n_words: int) -> np.ndarray:
        """Read n_words from consecutive records starting at start_record_index.

        Reads the 6 ENDF numeric fields from each record as a flat stream.
        """

        if n_words <= 0:
            return np.zeros((0,), dtype=float)

        out: list[float] = []
        i = start_record_index
        while len(out) < n_words:
            if i >= len(self.records):
                break
            rec = self.records[i]
            for x in rec.w:
                if len(out) >= n_words:
                    break
                out.append(float(x or 0.0))
            i += 1

        return np.asarray(out, dtype=float)

    def extract_mf3_mt(self, mt: int, itemp: int, nl: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """Extract MF=3 reaction data and scalar flux for a given temperature block.

        itemp and nl are 0-based indices.

        Returns (xs, flux) arrays shaped (ng, nsig0). For reactions that are missing,
        returns zeros.
        """

        ng = self.get_ng()
        nsig0 = self.get_nsig0()
        xs = np.zeros((ng, nsig0), dtype=float)
        flux = np.zeros((ng, nsig0), dtype=float)

        starts = self._section_starts(mf=3, mt=mt)
        if itemp < 0 or itemp >= len(starts):
            return xs, flux

        # NJOY GENDF MF=3 sections follow the pattern used in the legacy Fortran:
        #   - section header record (ns==1) contains nlgn=L1 and nsig0=L2
        #   - then repeating blocks of:
        #       control record (contains group index in word 6)
        #       fixed-length data block of size (2*nlgn*nsig0) words
        hdr_idx = starts[itemp]
        hdr = self.records[hdr_idx]
        nlgn = max(1, hdr.l1)
        nsig0_local = max(1, hdr.l2)
        if nsig0_local != nsig0:
            # Some MF=3 sections may not be sigma0-dependent; handle by using the section's nsig0.
            nsig0 = nsig0_local
            xs = np.zeros((ng, nsig0), dtype=float)
            flux = np.zeros((ng, nsig0), dtype=float)

        if nl < 0 or nl >= nlgn:
            raise ValueError(f"Requested Legendre nl={nl} but valid range is [0, {nlgn - 1}]")

        i = hdr_idx + 2  # first data record (control record is i-1)
        expected = 2 * nlgn * nsig0
        while i < len(self.records):
            rec = self.records[i]
            if rec.mf != 3 or rec.mt != mt:
                break
            # new section begins
            if rec.ns == 1:
                break

            ctrl = self.records[i - 1]
            group = int(ctrl.n2)
            if group <= 0 or group > ng:
                raise ValueError(f"MF3 MT{mt} itemp={itemp}: invalid group index {group} at record {i-1}")

            data = self.read_words(i, expected)
            if len(data) < expected:
                raise ValueError(
                    f"MF3 MT{mt} itemp={itemp}: not enough words at record {i} "
                    f"(got {len(data)}, expected {expected})."
                )

            flux_part = data[: nlgn * nsig0]
            xs_part = data[nlgn * nsig0 : 2 * nlgn * nsig0]

            for isig0 in range(nsig0):
                flux[group - 1, isig0] = flux_part[isig0 * nlgn + nl]
                xs[group - 1, isig0] = xs_part[isig0 * nlgn + nl]

            i = i + math.ceil(expected / 6) + 1  # skip data lines + next control record

        return xs, flux

    def extract_mf5_mt18_chi(self) -> np.ndarray:
        """Extract chi spectrum (MF=5 MT=18). Returns zeros if missing."""

        ng = self.get_ng()
        starts = self._section_starts(mf=5, mt=18)
        if not starts:
            return np.zeros((ng,), dtype=float)

        # For fission spectrum, NJOY typically writes a single section.
        hdr_idx = starts[0]
        i = hdr_idx + 1
        if i >= len(self.records):
            return np.zeros((ng,), dtype=float)

        rec = self.records[i]
        if rec.mf != 5 or rec.mt != 18:
            return np.zeros((ng,), dtype=float)

        n_words = rec.n1
        if n_words <= 0:
            return np.zeros((ng,), dtype=float)

        data = self.read_words(i + 1, n_words)
        out = np.zeros((ng,), dtype=float)
        out[: min(ng, len(data))] = data[: min(ng, len(data))]
        return out

    def extract_mf6_scatter(
        self, mt: int, itemp: int, nl: int = 0, *, strict: bool = True
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract MF=6 scattering for a given mt and temperature index.

        itemp and nl are 0-based indices.

        Returns (f, t, sig) where:
        - f: from-group indices (1-based)
        - t: to-group indices (1-based)
        - sig: (nonz, nsig0) values for requested Legendre nl
        """

        ng = self.get_ng()
        nsig0_global = self.get_nsig0()

        starts = self._section_starts(mf=6, mt=mt)
        if itemp < 0 or itemp >= len(starts):
            return (
                np.zeros((0,), dtype=int),
                np.zeros((0,), dtype=int),
                np.zeros((0, nsig0_global), dtype=float),
            )

        sec_idx = starts[itemp]
        first = self.records[sec_idx]

        # Determine nlgn and section-local nsig0:
        # - If the first record is a section header (N1==0), nlgn=L1 and nsig0=L2.
        # - Otherwise the first record is a control record; infer nlgn from N1 = nlgn*nsig0_global*ng2.
        if first.n1 == 0:
            nlgn = max(1, first.l1)
            nsig0_section = max(1, first.l2)
            i = sec_idx + 1
        else:
            # No explicit header record. Try to infer nlgn assuming this is a control record.
            ng2 = max(1, first.l1)
            n_words = first.n1
            nsig0_section = nsig0_global
            denom = nsig0_section * ng2
            if denom <= 0 or n_words % denom != 0:
                if strict:
                    raise ValueError(
                        f"MF6 MT{mt} itemp={itemp}: cannot infer nlgn from N1={n_words}, nsig0={nsig0_section}, ng2={ng2}"
                    )
                return (
                    np.zeros((0,), dtype=int),
                    np.zeros((0,), dtype=int),
                    np.zeros((0, nsig0_global), dtype=float),
                )
            nlgn = max(1, n_words // denom)
            i = sec_idx

        if nl < 0 or nl >= nlgn:
            raise ValueError(f"Requested Legendre nl={nl} but valid range is [0, {nlgn - 1}]")

        f_list: list[int] = []
        t_list: list[int] = []
        sig_list: list[np.ndarray] = []

        while i < len(self.records):
            rec = self.records[i]
            if rec.mf != 6 or rec.mt != mt:
                break
            # next temperature section begins
            if rec.ns == 1 and i != sec_idx:
                break

            ng2 = rec.l1
            ig2lo = rec.l2
            n_words = rec.n1
            ig = rec.n2

            if n_words <= 0 or ng2 <= 0:
                i += 1
                continue

            expected = nlgn * nsig0_section * ng2
            if n_words != expected:
                if strict:
                    raise ValueError(
                        f"MF6 MT{mt} itemp={itemp}: unexpected N1={n_words} (expected {expected}) at record {i}"
                    )
                return (
                    np.zeros((0,), dtype=int),
                    np.zeros((0,), dtype=int),
                    np.zeros((0, nsig0_global), dtype=float),
                )

            data = self.read_words(i + 1, n_words)
            if len(data) < expected:
                raise ValueError(
                    f"MF6 MT{mt} itemp={itemp}: not enough words at record {i} "
                    f"(got {len(data)}, expected {expected})."
                )

            # Skip flux block for this control record
            data = data[nlgn * nsig0_section :]

            k = 0
            # NJOY: ng2 includes the diagonal/self slot; Fortran loops to ng2-1
            for ito in range(ig2lo, ig2lo + ng2 - 1):
                f_list.append(ig)
                t_list.append(ito)

                row_section = np.zeros((nsig0_section,), dtype=float)
                for isig0 in range(nsig0_section):
                    row_section[isig0] = data[k + isig0 * nlgn + nl]
                if nsig0_section == nsig0_global:
                    row = row_section
                elif nsig0_section == 1:
                    row = np.full((nsig0_global,), row_section[0], dtype=float)
                else:
                    raise ValueError(
                        f"MF6 MT{mt} itemp={itemp}: cannot upcast section nsig0={nsig0_section} to global nsig0={nsig0_global}"
                    )
                sig_list.append(row)
                k += nlgn * nsig0_section

            i += 1 + math.ceil(n_words / 6)

        f_arr = np.asarray(f_list, dtype=int)
        t_arr = np.asarray(t_list, dtype=int)
        sig_arr = np.vstack(sig_list) if sig_list else np.zeros((0, nsig0_global), dtype=float)
        return f_arr, t_arr, sig_arr

    def _infer_mf6_nlgn(self, mt: int) -> int:
        # Look for a section header record where N1==0.
        for i in self._section_starts(mf=6, mt=mt):
            rec = self.records[i]
            if rec.n1 == 0:
                return max(1, rec.l1)
        return 1
