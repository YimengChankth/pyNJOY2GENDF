import openmc
import os
import glob
import numpy as np
import csv

def macroxsfiles2openmc(dir:str='material_xs', max_pn:int=0, export_to_h5:bool=True) -> openmc.MGXSLibrary :
    '''
    dir:str
        Folder where cross-section .csv files are located

    max_pn:int
        Number of legendre expansions to use

    export_to_h5:bool
        If True, then after generating library, save to {dir}/mgxs.h5

    '''

    # get the number of energy groups (the specific boundaries are unimportant)

    mat_list = get_material_names(dir)

    print(f'List of materials: {mat_list}')

    # read num energy groups
    with open(f'{dir}/{mat_list[0]}/macro/total','r') as f:
        ng = len(f.readlines()[0].split())
        groups = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, ng+1))

    # create library 
    mg_cross_sections_file = openmc.MGXSLibrary(groups)

    rxns1D = ['total','nu_fission', 'fission', 'chi', 'absorption',]

    for matname in mat_list:
        xsdata = openmc.XSdata(matname, groups)
        xsdata.order=max_pn

        for rxn in rxns1D:
            p = f'{dir}/{matname}/macro/{rxn}'
            setattr(xsdata, f'set_{rxn}', np.genfromtxt(p, dtype=float,delimiter=',', skip_header=0))

        scatter_xs = convert_sparse_csv_to_matrix(f'{dir}/{matname}/macro/scatter')
        # convert from LEGENDRE, FROM, TO to openmc's FROM, TO, LEGENDRE
        scatter_xs = scatter_xs.transpose(1, 2, 0)

        # limit to max pn order
        scatter_xs = scatter_xs[:,:,0:max_pn+1]

        # old scattering implementation, changed in favour of sparse representation
        # scatter_files=sorted(glob.glob(f'{dir}/{matname}/scatter_*'))

        # Pn = len(scatter_files)-1
        
        # # if there is not enough scatter files, then set max_pn to how many there are
        # max_pn = min([max_pn, Pn])

        # scatter_xs = np.zeros((ng, ng, max_pn+1))

        # for inl in range(max_pn+1):
        #     with open(f'{dir}/{mat_list[0]}/scatter_{inl}','r') as f:
        #         lines = f.readlines()
        #         for gfrom, line in enumerate(lines):
        #             scatter_xs[gfrom,:,inl] = [float(i) for i in line.split(',')]

        xsdata.set_scatter_matrix(scatter=scatter_xs)
        # Add the xs data to library
        mg_cross_sections_file.add_xsdata(xsdata)

    # write to disk
    if export_to_h5:
        mg_cross_sections_file.export_to_hdf5(f'{dir}/mgxs.h5')
        print(f'Saved mgxs file to {dir}/mgxs.h5')

    mg_cross_sections_file

def get_material_names(dir):
    assert os.path.exists(f'{dir}/material_list')
    with open(f'{dir}/material_list','r') as f:
        mat_list = f.readlines()[1:]
        mat_list = [i.strip() for i in mat_list]
    return mat_list

def convert_sparse_csv_to_matrix(sparse_csv:str):
    '''Convert a csv to numpy nd.array. The last column should be the value 
    '''
    rows = []
    with open(sparse_csv, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)

        for line in reader:
            # Last column is value
            val = float(line[-1])

            # All other columns are indices
            idx = tuple(int(x) for x in line[:-1])
            rows.append((idx, val))

    # Determine dimension sizes
    dims = []
    num_index_cols = len(rows[0][0])
    for col in range(num_index_cols):
        max_index = max(r[0][col] for r in rows)
        dims.append(max_index)

    # Create full array
    A = np.zeros(tuple(dims), dtype=float)

    # Fill array, convert from 1 based to 0 based indexing
    for idx, val in rows:
        zero_idx = tuple(i - 1 for i in idx)
        A[zero_idx] = val

    return A

def read_microxs(dir='material_xs') -> dict:
    '''Read microscopic cross-sections for all materials and return as dict
    '''

    mat_list = get_material_names(dir)

    microxsLib = {}

    for matname in mat_list:
        
        microxsLib[matname] = {}
        
        isolist = np.genfromtxt(fname=f'{dir}/{matname}/isolist', delimiter=',',skip_header=1, dtype=('U10', float, float))

        microxsLib[matname]['isolist'] = isolist

        xs1D = ['fission', 'chi', 'kerma', 'nu']
        xs2D = ['scatter_elas', 'scatter_inelas']

        for i in range(len(isolist)):
            isoname = str(isolist[i][0]).strip()
            microxsLib[matname][isoname] = {}

            for xs in xs1D:
                p =  f'{dir}/{matname}/{isoname}/{xs}'
                microxsLib[matname][isoname][xs] = np.genfromtxt(p, delimiter=',', skip_header=0)
            
            for xs in xs2D:
                p =  f'{dir}/{matname}/{isoname}/{xs}'
                microxsLib[matname][isoname][xs] = np.genfromtxt(p, delimiter=',',skip_header=1)

    microxsLib


if __name__ == "__main__":
    # macroxsfiles2openmc('new_ver_example/material_xs')
    read_microxs('new_ver_example/material_xs')
    pass