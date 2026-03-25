import numpy as np
import warnings

def print_tuples(tuples_list, break_every=3, rjust_size=12, lineprefix='   '):
    output = lineprefix
    for i, (x, y) in enumerate(tuples_list):
        output += f'{x}'.rjust(rjust_size) + f'{y}'.rjust(rjust_size)

        if i == (len(tuples_list) - 1):
            break  # avoid adding extra line break at the end

        if (i + 1) % break_every == 0:
            output += f'\n{lineprefix}'  # Add a line break after every 3 tuples

    return output

def read_endf_mt_mf(
    filename: str,
    MT: int,
    MF: int
) -> np.ndarray:
    """
    Read an ENDF file and extract all lines matching given MT and MF.

    Parameters
    ----------
    filename : str
        Path to ENDF file
    MT : int
        Reaction number
    MF : int
        File number

    Returns
    -------
    list[list[float]] or None
        Parsed data rows (6 floats per row) or None if not found
    """

    matched_rows = []

    with open(filename, "r") as f:
        for line in f:
            if len(line) < 75:
                continue  # skip malformed lines

            try:
                mf = int(line[70:72])
                mt = int(line[72:75])
            except ValueError:
                continue  # non-data line

            if mf == MF and mt == MT:
                values = []
                for i in range(6):
                    field = line[i * 11:(i + 1) * 11].strip()
                    if field:
                        # ENDF ends with +0, with no 'E', character. Manually do it ourself by finding either + or - character
                        idx = next((i for i, c in enumerate(field[::-1]) if c in ['+', '-']), -1)
                        if idx == -1:
                            # no +/- sign found
                            num = 0.0
                        else:
                            # split by number and exponent, but remember to reverse the index
                            idx = len(field) - idx - 1
                            num = float(f'{field[:idx]}E{field[idx:]}')  # test conversion
                        values.append(num)
                    else:
                        values.append(0.0)

                matched_rows.append(values)

    if not matched_rows:
        warnings.warn(
            f"MT={MT}, MF={MF} not found in ENDF file: {filename}",
            RuntimeWarning,
        )
        return None

    return np.array(matched_rows)