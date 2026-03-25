import numpy as np
import matplotlib.pyplot as plt
import openmc
import pickle
import matplotlib.pyplot as plt
import os
import csv

class kinfSolver:
    '''Solve the k-inf problem with Wielandt shift acceleration
    '''
    def __init__(self, total, scatter, nufission):

        total = np.asarray(total)
        scatter = np.asarray(scatter)
        nufission = np.asarray(nufission)


        self.total = np.diag(total.flatten())
        self.scatter = scatter # (Egroup from, Egroup to) matrix
        self.nufission = nufission # (Egroup from, Egroup to) matrix
        self.Egroups = len(total)

    @staticmethod
    def nufissionFromFissionChiNu(fission, chi, nu=None):
        '''Formulate the nu-fission matrix from fission, chi and nu

        nu-fission matrix are indexed in [Egroup from, Egroup to]

        '''

        fission = np.asarray(fission)
        fission = fission.flatten()
        chi = np.asarray(chi)
        chi = chi.flatten()
        if nu is None:
            nu = np.ones(len(chi))
        
        nu = np.asarray(nu)
        nu = nu.flatten()
        
        nufission = (nu*fission)[:,np.newaxis]*chi[np.newaxis,:]
        return nufission

    def solve(self, tol = 1E-6, verbose=True):
        '''Solve the k-inf problem using Wielandt shift (see: Space-Dependent Wielandt shifts for Multigroup Diffusion Eigenvalue Problems, Ben et al. 2017, Nuclear Science and Engineering)

        '''

        residual = 1
        iter = 0

        flux = np.ones((self.Egroups,1))/self.Egroups

        klambda = 1 # 1/kinf

        klambda_shift = 0 # Wielandt shift value

        # Wielandt shift parameters see: Wielandt Shifts for multigroup diffusion eigenvalue problems Yee et al.
        c1 = 10
        c0 = 0.01
        klambda_min = 1/3 #1/3

        max_iter = 1E5
        klambda_shift = 0.0

        while residual > tol:

            # Calculate fission + scattering source

            RHS = (klambda - klambda_shift)*np.matmul(self.nufission.T, flux)

            LHS = self.total - self.scatter.T - klambda_shift*self.nufission.T

            flux_new = np.linalg.solve(LHS, RHS)

            # Calculate new lambda
            temp = np.sum(np.matmul(self.nufission.T,klambda_shift*flux_new + (klambda - klambda_shift)*flux))
            klambda_new = temp/np.sum(np.matmul(self.nufission.T, flux_new))

            # Calculate new Wielandt shift, see "Space dependent
            # Wielandt Shifts for multigroup diffusion eigenvalue problems" Yee et
            # al. This is the implementation in PARCS            

            klambda_shift = max(klambda_new - c1*abs(klambda_new - klambda) - c0, klambda_min)

            # Calculate residual
            klambda_residual = abs(klambda_new - klambda)/klambda

            flux_residual = np.linalg.norm(np.abs(flux_new - flux)/flux)

            residual = max([klambda_residual, flux_residual])

            if verbose is True:
                print(f'Iteration: {iter}')
                print(f'\tkinf: {(1/klambda):1.3e}')
                print(f'\tResidual: {residual:1.3e}\n')

            # normalize fluxes 
            flux = flux_new/np.sum(flux_new)
            klambda = klambda_new

            iter = iter + 1
            if iter == max_iter:
                raise ValueError(f'Iteration number is {iter} and is outside of maximum iteration {max_iter}')

        if verbose is True:
            print(f'kinf calculation completed in {iter} iterations. \n\t Final kinf value: {1/klambda:1.6e} \n\t Final residual: {residual:1.6e} \n\n')
        
        self.flux = flux
        self.kinf = 1/klambda

        return flux.flatten(), 1/klambda

    @staticmethod
    def collapseFluxWeightedCrossSections(flux, XS, EBroadGroups):
        '''Collapse cross-sections using a specified flux and Energy group edges

        Parameters
        ----------
        flux : np.ndarray 
            The fine group flux distribution

        XS : 1D (total) or 2D np.ndarray (scatter, chinufission)
            XS to be collapsed

        EBroadGroups : np.ndarray
            The number of fine energy groups in each broad group. npte: sum(EBroadGroups) == number of energy groups

        '''

        numFineGroups = len(flux)
        if numFineGroups != np.sum(EBroadGroups):
            raise ValueError(f'Number of groups in flux is: {numFineGroups}, but sum(EbroadGroups) = {np.sum(EBroadGroups)}')

        numBroadGroups = len(EBroadGroups)

        EBroadGroups = np.cumsum(EBroadGroups)

        # Calculate the broadGroupFluxes

        broadGroupFlux = np.zeros(numBroadGroups)

        for g in range(numBroadGroups):
            if g==0:
                gstart = 0
            else:
                gstart = EBroadGroups[g-1]

            broadGroupFlux[g] = np.sum(flux[gstart:EBroadGroups[g]])

        # Check if XS is 1D or 2D
        
        if len(XS.shape) == 1:
            collapsedXS = np.zeros(numBroadGroups)
            # XS is 1D, e.g. total cross-section
            for g in range(numBroadGroups):
                if g==0:
                    gstart = 0
                else:
                    gstart = EBroadGroups[g-1]
                rxnRate = np.dot(flux[gstart:EBroadGroups[g]], XS[gstart:EBroadGroups[g]])
                collapsedXS[g] = rxnRate/broadGroupFlux[g]

            return collapsedXS


        if len(XS.shape) == 2:
            # XS is 2D e.g. scatter matrix, chi-nu-fission matrix. Format:[Gfrom, Gto]

            # Calculate sum of reaction rates
            temp = np.zeros((numFineGroups, numBroadGroups))
            for gfrom in range(numFineGroups):
                for gto in range(numBroadGroups):
                    if gto == 0:
                        gstart = 0
                    else:
                        gstart = EBroadGroups[gto-1]
                    gend = EBroadGroups[gto]

                    temp[gfrom, gto] = np.sum(XS[gfrom, gstart:gend])*flux[gfrom]


            # Flux weighting
            collapsedXS = np.zeros((numBroadGroups, numBroadGroups))
            for gfrom in range(numBroadGroups):
                if gfrom == 0:
                    gstart = 0
                else:
                    gstart = EBroadGroups[gto-1]

                gend = EBroadGroups[gfrom]

                collapsedXS[gfrom,:] = np.sum(temp[gstart:gend,:], axis=0)/broadGroupFlux[gfrom]
            
            return collapsedXS

        pass

    @staticmethod
    def collapseChi(chi, EBroadGroups):
        '''Collapse Chi. There is no flux weighting 

        Parameters
        ----------
        chi : 1D (total) or 2D np.ndarray (scatter, chinufission)
            XS to be collapsed

        EBroadGroups : np.ndarray
            The number of fine energy groups in each broad group. npte: sum(EBroadGroups) == number of energy groups

        '''
        numFineGroups = len(chi)
        if numFineGroups != np.sum(EBroadGroups):
            raise ValueError(f'Number of groups in flux is: {numFineGroups}, but sum(EbroadGroups) = {np.sum(EBroadGroups)}')

        colChi = np.zeros(len(EBroadGroups))
        count = 0
        for g in range(len(EBroadGroups)):
            colChi[g] = np.sum(chi[count:count+EBroadGroups[g]])
            count = count + EBroadGroups[g]

        return colChi

def run_kinfsolver(material:openmc.XSdata, temperature_index:int=0):
    '''
    Runs infinite medium problem for given openmc XSdata object. Writes to a results file the flux, keff, absorption rate, fission rate 
    
    :param material: Description
    :type material: openmc.XSdata

    temperature_index : int
        which temperature index to use. As a material can have temperature dependent XSs.

    '''

    total = material.total[temperature_index]
    scatter = material.scatter_matrix[temperature_index][:,:,0]
    chi = material.chi[temperature_index]
    nu_fission = material.nu_fission[temperature_index]

    chinufiss = kinfSolver.nufissionFromFissionChiNu(fission=nu_fission, chi=chi,nu=None)

    kinfsolver = kinfSolver(
                            total=total,
                            scatter = scatter,
                            nufission = chinufiss
                            )
    

    flux, kinf = kinfsolver.solve()

    # perform flux normalization 
    flux = flux/np.sum(flux)

    absorption = total - np.sum(scatter, axis=1)

    abs_rate = np.multiply(absorption, flux)
    fisssource_rate = np.multiply(nu_fission, flux)

    results = {
        'flux':flux, 
        'keff':kinf,
        'abs_rate':abs_rate,
        'fisssource_rate':fisssource_rate,
        'macro_nufiss':nu_fission,
        'macro_abs':absorption

    }

    return results

def plot_results(results, outdir):
    '''
    Docstring for plot_results
    
    :param results: Description
    '''
    os.makedirs(outdir, exist_ok=True)
    
    def plot_vec(vals, title:str, ylabel:str, savepath:str):
        fig = plt.Figure()
        plt.plot(vals)
        plt.yscale('log')
        plt.ylabel(ylabel)
        plt.xlabel('Energy Group index')
        plt.title(title)
        plt.tight_layout()
        plt.savefig(savepath)
        plt.close()

    keff = results['keff']
    title=f'keff: {keff:1.6f}'
    plot_vec(results['flux'], title, r'$\phi_{g} [1/cm^{2}-s]$', f'{outdir}/flux.pdf')
    plot_vec(results['fisssource_rate'], title, r'$\nu_{g}\Sigma_{f,g}\phi_{g} [1/cm^{3}-s]$', f'{outdir}/fisssource_rate.pdf')
    plot_vec(results['abs_rate'], title, r'$\Sigma_{a,g}\phi_{g} [1/cm^{3}-s]$', f'{outdir}/abs_rate.pdf')

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
    # data = np.genfromtxt(f'{dir}/{mat_list[0]}/macro/total.csv', delimiter=',', skip_header=1)
    ng = len(convert_csv_to_matrix(f'{dir}/{mat_list[0]}/macro/total.csv'))

    # This is arbitrary
    groups = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, ng+1))

    # create library 
    mg_cross_sections_file = openmc.MGXSLibrary(groups)

    rxns1D = ['total','nu_fission', 'fission', 'chi', 'absorption',]

    for matname in mat_list:
        xsdata = openmc.XSdata(matname, groups)
        xsdata.order=max_pn

        for rxn in rxns1D:
            p = f'{dir}/{matname}/macro/{rxn}.csv'
            vals = convert_csv_to_matrix(p)
            getattr(xsdata, f'set_{rxn}')(vals) # np.genfromtxt(p, dtype=float,delimiter=',', skip_header=0)
            
        scatter_xs = convert_csv_to_matrix(f'{dir}/{matname}/macro/scatter.csv')
        # convert from LEGENDRE, FROM, TO to openmc's FROM, TO, LEGENDRE
        scatter_xs = scatter_xs.transpose(1, 2, 0)

        # limit to max pn order
        scatter_xs = scatter_xs[:,:,0:max_pn+1]

        xsdata.set_scatter_matrix(scatter=scatter_xs)
        # Add the xs data to library
        mg_cross_sections_file.add_xsdata(xsdata)

    # write to disk
    if export_to_h5:
        mg_cross_sections_file.export_to_hdf5(f'{dir}/mgxs.h5')
        print(f'Saved mgxs file to {dir}/mgxs.h5')

    return mg_cross_sections_file

def get_material_names(dir):
    assert os.path.exists(f'{dir}/material_list')
    with open(f'{dir}/material_list','r') as f:
        mat_list = f.readlines()[1:]
        mat_list = [i.strip() for i in mat_list]
    return mat_list

def convert_csv_to_matrix(csv_path:str):
    '''Convert a csv to numpy nd.array. The last column should be the value 
    '''
    rows = []
    with open(csv_path, newline="") as f:
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


if __name__ == "__main__":

    # get mgxs object
    materials = macroxsfiles2openmc(dir='material_xs', export_to_h5=False)
    inner_material = materials.get_by_name('INNER')
    results = run_kinfsolver(inner_material, 0)
    with open(f'kinf_solve_results.pkl','wb') as f:
        pickle.dump(results, f)
    
    plot_results(results, outdir='result_plots')

    pass
