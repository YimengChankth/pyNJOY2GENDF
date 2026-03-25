# read wims_njoy_stc.txt

def read_purr_file():
    with open('wims_njoy_stc.txt', 'r') as f:
        lines = f.readlines()

    nuc_list = []

    for c,v in enumerate(lines):
        if (v.startswith('purr')):
            cline = c - 15
            if 'ENDF/B-VII' in lines[cline]:
                # print(f'cl: {cline}, line: {lines[cline]}')
                nuclide = lines[cline].split()[0]
                tmp = nuclide.split('-')
                if len(tmp) < 3:
                    continue
                iso = tmp[1]
                m = tmp[2]

                iso_str = f'{iso}{m}'
                iso_str = iso_str.capitalize()
                if iso_str not in nuc_list:
                    nuc_list.append(iso_str)
                    print(f"'{iso_str}',")

if __name__ == "__main__":
    read_purr_file()