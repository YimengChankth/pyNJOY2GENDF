import os 
import glob
import pickle
import re
import numpy as np
from pyNJOY2GENDF.card8 import WeightingFlux, GROUPRCard8
from pyNJOY2GENDF.templates import default_BROADR, default_GROUPR, default_HEATR, default_RECONR, default_UNRESR, default_PURR, RXN, GROUPRCard8
import shutil
import subprocess
from pyNJOY2GENDF.njoy_err_analysis import write_only_errors_to_file
import datetime
import glob
import warnings
from pyNJOY2GENDF.wims_params import WIMS_EHI_SIGPOT, WIMS_PURR

# parallel run of NJOY
from multiprocessing import Pool

class njoy2gendf:
    def __init__(self, savedir:str, path2endf:str, njoy_exec:str='njoy21', nuclide_list:list[str]=['U235', 'U238'],filename_convention='ENDF-B'):
        '''
        Docstring for __init__
        
        :param savedir: str
            Where to save the results 
        :param path2endf: str
            path to the folder containing endf files 
        :param njoy_exec: str
            Path to NJOY executable
        :param filename_convention : str 
            Either 'ENDF-B' style e.g. n-092_U_235.endf OR 'JEFF' style e.g. 92-U-235g.jeff33
        '''

        os.makedirs(savedir, exist_ok=True)
        self.savedir = savedir

        assert os.path.exists(path2endf), f'Path to ENDF files does not exist: {path2endf}'
        self.path2endf = path2endf
        self.njoy_exec = njoy_exec
        self.nuclide_list = nuclide_list

        assert filename_convention in ['ENDF-B', 'JEFF']
        self.filename_convention = filename_convention
        pass

    @staticmethod
    def loadisotope_matd(file='iso2matd_table.pkl'):
        assert os.path.exists(file)
        with open(file, 'rb') as f:
            iso2matd = pickle.load(f)

        return iso2matd

    @staticmethod
    def getisotope_matd(endf_files:list[str]) -> dict:
        result = []
        for f in endf_files:
            with open(f, 'r') as fp:
                line2 = fp.readlines()[1]
                matd = int(line2[66:70])
            result.append(matd)
        return result

    def write_njoy_gendf_inputs(self, 
                                t300_rxns:dict[str, list[RXN]]=None,
                                tother_rxns:dict[str, list[RXN]]=None,
                                energybin_edges=[1e-5, 0.025, 1e6],
                                library_name='ENDF-B/VIII',
                                temperatures = [300.,    600.,    900],
                                sig0s = [1.e10,   1.e5,    1.e4,    1000.,   100.,    10.,     1.],
                                use_PURR:list[str]=[], # if nuclide is in use_PURR then PURR: unresolved resonance data processing for isotopes will be used instead of UNRESR (default)
                                use_HEATR:bool=True,
                                card8_params:dict[str, GROUPRCard8] | None = None,
                                ):
        '''This writes NJOY inputs similar to what KM provided. 

        Modifications to this class (or inheritance) should overwrite the methods inbuilt 

        Order of inputs:
        RECONR
        BROADR
        UNRESR
        HEATR
        GROUPR

        Parameters
        ----------
        nuclide_list : list of str
            nuclide names

        t300_rxns : dict of nuclide to list of RXN objects for base temperature

        tother_rxns : dict of nuclide to list of RXN objects for other temperatures

        energy_bin_edges : list of float
            eV for default NJOY. note that SERPENT is default MeV

        library_name : str
            cosmetic, just to tell users where library is from by reading input file 

        temperatures : list of floats
            temperature to evaluate BROADR

        sig0s : list of floats
            background XS 

        card8a_params : dict of str to GROUPRCard8
            If a nuclide is in card8a_params, the corresponding GROUPRCard8 object will be used for that nuclide.


        '''


        # valdiate reactions
        if t300_rxns is None:
            print('No rxns provided for base temperature (t300_rxns). Using default reactions based on ENDF file MT numbers and includes basic reactions', flush=True)
            print(f'If you want to specify reactions, provide t300_rxns as a dict of nuclide to list of RXN objects for base temperature reactions. See get_T300_rxn method for details', flush=True)
            t300_rxns = {nuclide: self.get_T300_rxn(self.find_endf_file_for_nuclide(nuclide), use_HEATR=use_HEATR) for nuclide in self.nuclide_list}
        
        if ((tother_rxns is None) and (len(temperatures) > 1)):
            print('No rxns provided for non-base temperature (tother_rxns). Using default reactions based on ENDF file MT numbers and includes basic reactions. For higher temperatures, only total, elastic, fission, (n,gamma), nubar, and heat (if HEATR is used) are included based on NJOY manual recommendation', flush=True)
            print(f'If you want to specify reactions, provide tother_rxns as a dict of nuclide to list of RXN objects for other temperature reactions. See get_Tother_rxn method for details', flush=True)
            tother_rxns = {nuclide: self.get_Tother_rxn(self.find_endf_file_for_nuclide(nuclide)) for nuclide in self.nuclide_list}

        # Validate sig0s requirements
        if len(sig0s) == 0:
            raise ValueError("sig0s cannot be empty")
        
        # Check that sig0s are in strict descending order
        for i in range(len(sig0s) - 1):
            if sig0s[i] <= sig0s[i + 1]:
                raise ValueError(
                    f"sig0s must be in strict descending order. "
                    f"Found sig0s[{i}]={sig0s[i]} <= sig0s[{i+1}]={sig0s[i+1]}"
                )
        
        # Check that first sig0 value is 1.e10 (infinite dilution)
        if abs(sig0s[0] - 1.e10) > 1.e-6:
            raise ValueError(
                f"sig0s[0] must be 1.e10 for infinite dilution. "
                f"Got sig0s[0]={sig0s[0]}"
            )
        
        # if card8a_params is None, then use WIMS parameters 
        # if card8a_params is None:
        #     print('No card8a_params provided. Using WIMS parameters for EHI and SIGPOT for nuclides that are available in the WIMS library. For nuclides not in the WIMS library, default iwt will be used and card 8a will be skipped.', flush=True)
        #     card8a_params = {}
        #     for nuclide in self.nuclide_list:
        #         if nuclide in WIMS_EHI_SIGPOT.keys():
        #             card8a_params[nuclide] = GROUPRCard8()
        #             card8a_params[nuclide].fehi = WIMS_EHI_SIGPOT[nuclide][0]
        #             card8a_params[nuclide].sigpot = WIMS_EHI_SIGPOT[nuclide][1]
        #             card8a_params[nuclide].nflmax = 20000
        #             card8a_params[nuclide].weightingflux = groupr_weighting_flux
        #             card8a_params[nuclide].iwt = -1  # use -1 for flux weight calculator
        #         else:
        #             card8a_params[nuclide] = None

        # if use_PURR is None:
        #     use_PURR = WIMS_PURR


        os.makedirs(f'{self.savedir}/NJOY_INPUTS', exist_ok=True)

        # user input check for card 8
        # if groupr_iwt == 'weighting_flux':
        #     assert isinstance(groupr_weighting_flux, WeightingFlux), 'If groupr_iwt is "weighting_flux", groupr_weighting_flux must be provided as WeightingFlux object'
        #     assert groupr_fehi_sigpot_dict is not None, 'If groupr_iwt is "weighting_flux", groupr_fehi_sigpot_dict must be provided'
        # else:
        #     assert np.abs(groupr_iwt) in range(2,13), '|groupr_iwt| must be integer between 2 and 8" if not using weighting flux'
        #     assert groupr_weighting_flux is None, 'If |groupr_iwt| is 2-13, groupr_weighting_flux must be None'

        # # if no groupr_fehi then instantiate empty dict
        # if groupr_fehi_sigpot_dict is None:
        #     groupr_fehi_sigpot_dict = {}


        # iso2matd = self.loadisotope_matd('iso2matd_table.pkl')
        print(f'NJOY2GENDF info: ')
        print(f'Writing NJOY input to: {self.savedir}/NJOY_INPUTS. \n ENDF file read from: {self.path2endf}', flush=True)
        print(f'The following nucldies will be processed: {self.nuclide_list} \n', flush=True)

        purr_nucs = [nuclide for nuclide in self.nuclide_list if nuclide in use_PURR]

        print(f'The following nuclides will use PURR instead of UNRESR for unresolved resonance processing: {purr_nucs} \n', flush=True)

        for nuclide in self.nuclide_list:
            # Default tapes
            endf_tape = 20
            final_tape = 21

            print('='*50, flush=True)
            print(f'Writing input for nuclide: {nuclide}', flush=True)
            
            # Get matds from the libary (note that between ENDF and JEFF there are small differences)
            endf_file = self.find_endf_file_for_nuclide(nuclide)

            # if unable to find the file then raise error 
            assert endf_file is not None

            matd = self.getisotope_matd([endf_file])[0]
            header_comment = f'''-- Generated by njoy2GENDF.py on {datetime.datetime.now().strftime("%Y/%m/%d")} 
-- for {nuclide} with library located at {self.path2endf} \n\n'''
            
            reconr = self.generate_RECONR(nendf=endf_tape, npend=final_tape, nuclide=nuclide, mat=matd, library=library_name)

            broadr = self.generate_BROADR(nendf=endf_tape, nin=final_tape, nout=final_tape+1, mat=matd, temperatures=temperatures)
            final_tape += 1


            if nuclide in use_PURR:
                print(f'Using PURR unresolved resonance processing for nuclide {nuclide}', flush=True)
                unresr = self.generate_PURR(nendf=endf_tape, nin=final_tape, nout=final_tape+1, mat=matd, temperatures=temperatures, sig0s=sig0s)
                final_tape += 1
            else:
                print(f'Using UNRESR unresolved resonance processing for nuclide {nuclide}', flush=True)
                unresr = self.generate_UNRESR(nendf=endf_tape, nin=final_tape, nout=final_tape+1, mat=matd, temperatures=temperatures, sig0s=sig0s)
                final_tape += 1

            if use_HEATR:
                heatr = self.generate_HEATR(nendf=endf_tape, nin=final_tape, nout=final_tape+1, mat=matd)
                final_tape += 1
            else:
                heatr = ''

            # # get room temp reactions
            # t300_rxns = self.get_T300_rxn(nuclide, use_HEATR=use_HEATR)

            # # get temperature dependent reactions. See NJOY manual page 244.
            # tother_rxns = self.get_Tother_rxn(nuclide)

            card8_nuc = card8_params[nuclide]

            rxn_base = t300_rxns[nuclide]
            if tother_rxns is not None:
                rxn_temp = tother_rxns[nuclide]
            else:
                rxn_temp = None

            groupr = self.generate_GROUPR(
                                          nendf=endf_tape,
                                          npend=final_tape,
                                          ngout2=final_tape+1,
                                          nuclide=nuclide, 
                                          mat=matd, 
                                          library=library_name, 
                                          temperatures=temperatures, 
                                          sig0s=sig0s, 
                                          energybin_edges=energybin_edges, 
                                          rxn_base=rxn_base, 
                                          rxn_temp=rxn_temp, 
                                          groupr_iwt = card8_nuc.iwt,
                                          groupr_card8 = card8_nuc.write_to_block()
                                          )

            # addition of strings
            inp_file_lines = header_comment + reconr + broadr + unresr + heatr + groupr + 'STOP'

            # write to file!
            with open(f'{self.savedir}/NJOY_INPUTS/{nuclide}.inp','w') as f:
                f.write(inp_file_lines)

            print('\n\n', flush=True)

            pass
        pass
    
    def find_endf_file_for_nuclide(self,nuclide:str):
        '''
        Find the correct xs data file for given nuclide. This should be modified for JEFF etc. 
        
        :param self: Description
        :param nuclide: name of nuclide with convention e.g. Am241
        '''

        assert nuclide[0].isalpha(), "Must start with a letter"

        for c, v in enumerate(nuclide):
            if v.isdigit():
                el = nuclide[0:c]
                iso = nuclide[c:]
                break

        # ENDF-B style e.g. n-092_U_235.endf 
        found_flag = False
        if self.filename_convention == 'ENDF-B':
            za = f'_{el}_' + f'{iso}'.rjust(3, '0')
            endf_files = glob.glob(f'{self.path2endf}/*.endf')
            for file in endf_files:
                tmp = file.split('/')[-1][0:-5]
                if za in tmp:
                    found_flag = True
                    return file
        
        # JEFF naming convention style e.g. 92-U-235g. This works for JEFF4.0. 
        if self.filename_convention == 'JEFF':
            za = f'-{el}-' + f'{iso}'.rjust(3,'0')
            jeff_files = glob.glob(f'{self.path2endf}/*.jeff*')
            for file in jeff_files:
                tmp = file.split('/')[-1]
                
                if za in tmp:
                    found_flag = True
                    return file

        if found_flag == False:
            print(f'Nuclide {nuclide} not found in {self.path2endf}. If you expect the nuclide, modify the method `find_endf_file_for_nuclide()` in njoy2gendf.py', flush=True)

        pass
    
    @staticmethod
    def get_T300_rxn(endf_file, use_HEATR:bool=True):  
        '''
        Return a list of RXN objects for reactions for base temperature. 
    3       1       total / mfd mtd mtname
    3       102     ngamma /
    3       18      fission /
    8       16      n2n / 
    8       2       elastic /
    6       51      inelastic /
    6       -90     inelastic /

        :param endf_file: Path to the ENDF file
        :param use_HEATR: Whether to include HEATR reactions
        '''

        # check if nuclide exists in data file

        assert(endf_file is not None)

        # total reaction is always present 
        rxn_list = [RXN(mfd=3, mt=1, mtname='total')]
        
        # check (n,gamma)
        ng_mf = check_rxn_in_endf(mt=102, endf_file=endf_file)
        if (3 in ng_mf):
            rxn_list.append(RXN(mfd=3, mt=102, mtname='ngamma'))

        # check if fissile
        fiss_mf = check_rxn_in_endf(mt=18, endf_file=endf_file)

        if (3 in fiss_mf):
            # if fission present, add in other relevant reactions
            rxn_list += [RXN(mfd=3, mt=18, mtname='fission'),
                         RXN(mfd=3, mt=259, mtname='inv'),
                         RXN(mfd=3, mt=452, mtname='nubar'),
                         RXN(mfd=5, mt=18, mtname='chi'), # note that mfd=5 is the outgoing energy distribution. We assume isotropic fission. If no isotopic then request mfd=6 which also gets angular distribution out. Not relevant for fission neutron energies
                         ]

        # elastic

        elas_mf = check_rxn_in_endf(mt=2, endf_file=endf_file)
        if (6 in elas_mf):
            rxn_list += [RXN(mfd=8, mt=2, mtname='elastic')]
        else:
            # if mf=6 is not present then there must be mf=4 (4 is outgoing angular distribution)
            assert 4 in elas_mf
            rxn_list += [RXN(mfd=6, mt=2, mtname='elastic')]

        # (n2n) if present
        n2n_mf = check_rxn_in_endf(mt=16, endf_file=endf_file)
        if (6 in n2n_mf):
            rxn_list += [RXN(mfd=8, mt=16, mtname='n2n')]
        
        # inelastic. use notation e.g. 51, -71 to denote the processing of inelastic from 51 to 71, Will also output 91 as inelastic to continuum
        inelas = get_inelastic(endf_file)
        # if there is no inelastic then skip
        if inelas is not None:
            rxn_list += inelas

        if use_HEATR:
            rxn_list += [RXN(mfd=3, mt=301, mtname='heat')] # heating]

        return rxn_list

    @staticmethod
    def get_Tother_rxn(endf_file):
        '''
        Return a list of RXN objects for reactions for non T=300K temperatures. Accoring to NJOY manual page 249: 
        
        At higher temperatures, the threshold reactions should be
        omitted, because their cross sections do not change significantly with temperature
        except at the most extreme conditions. This means that only the following
        reactions should be included for the higher temperatures (if present): total
        (MT=1), elastic (MT=2), fission (MT=18), radiative capture (MT=102), heating
        (MT=301), kinematic KERMA (MT=443), damage (MT=444), and any
        thermal cross sections (MT=221-250). Only the elastic and thermal matrices
        should be included at the higher temperatures.

        For KM's MONTE code we get only the following MTs:
            3       1       total /
            3       102     ngamma /
            3       18      fission /
            6       2       elastic /
            3       452     nubar /


        :param nuclide:str
            Name of nuclide e.g. Am241
        '''

        # check if nuclide exists in data file
        assert(endf_file is not None)

        # total reaction is always present 
        rxn_list = [RXN(mfd=3, mt=1, mtname='total')]
        
        # check (n,gamma)
        ng_mf = check_rxn_in_endf(mt=102, endf_file=endf_file)
        if (3 in ng_mf):
            rxn_list.append(RXN(mfd=3, mt=102, mtname='ngamma'))

        # check if fissile
        fiss_mf = check_rxn_in_endf(mt=18, endf_file=endf_file)

        if (3 in fiss_mf):
            # if fission present, add in other relevant reactions
            rxn_list += [RXN(mfd=3, mt=18, mtname='fission'),
                         RXN(mfd=3, mt=452, mtname='nubar'),
                         ]

        # elastic - need MF=6 for scattering matrix at higher temperatures
        rxn_list += [RXN(mfd=6, mt=2, mtname='elastic')]

        return rxn_list

    def generate_RECONR(self,  **kwargs):
        reconr = default_RECONR.write_block(**kwargs)
        return reconr

    def generate_BROADR(self, **kwargs):
        broadr = default_BROADR.write_block(**kwargs)
        return broadr

    def generate_UNRESR(self, **kwargs):
        unresr = default_UNRESR.write_block(**kwargs)
        return unresr
    
    def generate_PURR(self, **kwargs):
        purr = default_PURR.write_block(**kwargs)
        return purr

    def generate_HEATR(self,  **kwargs):
        heatr = default_HEATR.write_block(**kwargs)
        return heatr

    def generate_GROUPR(self, **kwargs):
        groupr = default_GROUPR.write_block(**kwargs)
        return groupr

    def run_njoy_inputs(self, save_pendf=False, nthreads=1):
        '''
        Docstring for run_njoy_inputs. Run self.write_njoy_gendf_inputs() first!
        
        nthreads: number of threads to use in parallel NJOY runs.

        :param self: Description
        :param save_pendf: Whether to save the PENDF files
        '''

        assert os.path.exists(f'{self.savedir}/NJOY_INPUTS')

        inp_files = []
        # get all .inp files 
        for nuc in self.nuclide_list:
            p = f'{self.savedir}/NJOY_INPUTS/{nuc}.inp'
            if os.path.exists(p):
                inp_files.append(p)
            else:
                print(f'WARNING: nuclide {nuc} no input found in {self.savedir}/NJOY_INPUTS. Check that write_njoy_gendf_inputs method has been run', flush=True)

        # folder for njoy output
        out_dir = f'{self.savedir}/NJOY_OUTPUTS'
        os.makedirs(out_dir, exist_ok=True)

        # folder for GENDF output
        gendf_dir = f'{self.savedir}/GENDF'
        os.makedirs(gendf_dir, exist_ok=True)

        # folder for PENDF output if save_pendf is True
        if save_pendf:
            pendf_dir = f'{self.savedir}/PENDF'
            os.makedirs(pendf_dir, exist_ok=True)

        # it might make sense to run NJOY in a temp folder to avoid any clashes with tapeXX files 
        # perform paralel runs here
        endf_files = []
        for inp in inp_files:
            nuclide = inp.split('/')[-1][0:-4]
            # copy over ENDF file as tape20
            endf_files.append(self.find_endf_file_for_nuclide(nuclide))

        inp_endf_pairs = list(zip(inp_files, endf_files))

        # Create argument tuples for each call
        args_list = [(pair, save_pendf, self.savedir, out_dir, self.njoy_exec) for pair in inp_endf_pairs]

        with Pool(processes=nthreads) as pool:
            pool.starmap(self.run_single_njoy_input, args_list)

        self.check_njoy_runs()

    @staticmethod
    def run_single_njoy_input(inp_endf_pairs, save_pendf, base_savedir, outdir, njoy_exec='njoy21'):
        '''
        Docstring for run_single_njoy_input
        :param inp_endf_pairs: tuple of (inp_file, endf_file)
        :param save_pendf: Whether to save the PENDF files (tape 24)
        :param outdir: Description
        :param inp_file: Description
        :param njoy_exec: Description
        '''

        inp_file = inp_endf_pairs[0]
        
        print(f'Running NJOY for {inp_file}', flush=True)

        nuclide = inp_file.split('/')[-1][0:-4]

        os.makedirs(f'tmp_{nuclide}',exist_ok=True)

        # copy over ENDF file as tape20
        shutil.copy(inp_endf_pairs[1], f'tmp_{nuclide}/tape20')
        # copy over .inp file as njoy inpuy
        
        # clear old output (NJOY will append otherwise)
        if os.path.exists(f'{outdir}/{nuclide}.out') == True:
            os.remove(f'{outdir}/{nuclide}.out')

        # run njoy
        shutil.copy(inp_file, f'tmp_{nuclide}/njoy_input.inp')
        os.chdir(f'tmp_{nuclide}')
        # run njoy, suppress output
        subprocess.run(
            [
                njoy_exec,
                "-i", "njoy_input.inp",
                "-o", f"../{outdir}/{nuclide}.out", 
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            #check=True, # disable check so that even if njoy fails we can still capture the output files
        )
        # if successful, then several files will be created all beginning with tapeXX
        # copy the last tape with the largest integer to gendf folder   

        # find largest tape number in tapeXX
        tape_files = glob.glob('tape*')
        # Extract numbers from filenames
        numbers = []
        for file in tape_files:
            match = re.search(r'tape(\d+)', file)
            if match:
                numbers.append(int(match.group(1)))

        # Find the largest number,this will be the final GENDF file
        max_n = max(numbers)
        shutil.copy(f'tape{max_n}', f'../{base_savedir}/GENDF/{nuclide}.gendf')

        # the second last is the pendf file, save if save_pendf is True
        if save_pendf:
            if len(numbers) >= 2:
                second_last_n = sorted(numbers)[-2]
                shutil.copy(f'tape{second_last_n}', f'../{base_savedir}/PENDF/{nuclide}.pendf')
            else:
                print(f'WARNING: Only one tape file found for {nuclide}. Cannot save PENDF file.', flush=True)

        # clean up 
        tape_files = glob.glob('tape*')
        for tape in tape_files:
            if os.path.exists(tape):
                os.remove(tape)

        # card8 flux values (if present) (if ninwt is = -10 or -20)
        flux_files =  glob.glob('calculated_flux.*')
        os.makedirs(f'../{outdir}/NJOY_FLUXCALC', exist_ok=True)
        for ff in flux_files:
            if os.path.exists(ff):
                shutil.copy(ff, f'../{outdir}/NJOY_FLUXCALC/{nuclide}.{ff}')

        os.chdir('..')
        shutil.rmtree(f'tmp_{nuclide}')


    def check_njoy_runs(self):
        '''
        Produces a text file (njoy_run_analysis.tex, open as text file ) in latex table format showing if there are any others in the folder self.savedir/NJOY_OUTPUTS
        
        :param self: Description
        '''
        assert os.path.exists(f'{self.savedir}/NJOY_OUTPUTS')

        output_table = write_only_errors_to_file(f'{self.savedir}/NJOY_OUTPUTS', self.nuclide_list)

        with open(f'{self.savedir}/NJOY_OUTPUTS/njoy_run_analysis.tex','w') as f:
            f.write(output_table)

        print(f'Wrote NJOY analysis to: {self.savedir}/NJOY_OUTPUTS/njoy_run_analysis.tex')

def check_rxn_in_endf(mt:int, endf_file):
    '''for a given reaction, check if present in endf_file. If found, return also the mf number 
    '''

    mf_mt = parse_endf_mf_mt(endf_file)

    mf_number = []
    for mf_present in mf_mt.keys():
        if mt in mf_mt[mf_present]:
            mf_number.append(mf_present)
        
    return mf_number


def parse_endf_mf_mt(filename):
    mf_mt = {
        3:set(),
        4:set(),
        5:set(),
        6:set(),
    }  # {MF: set(MT numbers)}

    with open(filename, "r") as f:
        for line in f:
            if len(line) < 75:
                continue  # skip malformed lines

            try:
                mf = int(line[70:72])
                mt = int(line[72:75])
            except ValueError:
                continue  # skip lines that don't contain valid integers

            if mf not in mf_mt:
                mf_mt[mf] = set()
            mf_mt[mf].add(mt)

    for key, item in mf_mt.items():
        tmp = list(item)
        tmp.sort()
        mf_mt[key] = tmp


    return mf_mt    

def get_inelastic(endf_file) -> list['RXN', 'RXN'] | list['RXN', 'RXN', 'RXN']:

    '''
    Returns inelastic RXN objects, checking for MT numbers 51-91. 
    E.g. say an ENDF file has reaction number [51,52,...86, 91]
    6 51 inelastic
    6 -86 inelastic
    6 91 inelastic(continuum)

    :param endf_file: path to endf
    :return: 
    :rtype: tuple
    '''

    mf_mt = parse_endf_mf_mt(endf_file)

    # check if reaction 51 is in mf=4/5 or 6
    found_inel = False
    if (51 in mf_mt[4]) or (51 in mf_mt[5]):
        found_inel = True
        if 51 in mf_mt[4]:
            inel_mt = 4
        else:
            inel_mt = 5
    
    if 51 in (mf_mt[6]):
        found_inel = True
        inel_mt = 6
    
    # if no inel
    if found_inel == False:
        return None
    
    tmp = np.array(mf_mt[inel_mt])
    inel_list = tmp[(tmp >= 51) & (tmp <= 90)]
    min_inel = np.min(inel_list)
    rxn_min = RXN.from_mf(mf=inel_mt, mt=min_inel, mtname='inelastic')
    inel_list = list(inel_list)
    max_inel =  np.max(inel_list)
    rxn_max = RXN.from_mf(mf=inel_mt, mt=-max_inel, mtname='inelastic') # turn negative  (-N) so njoy will know to process from 51 to N
    
    #  check for continuum 
    if (91 in mf_mt[4]) or (91 in mf_mt[6]):
        # 91 found 
        if (91 in mf_mt[4]):
            inel_mt_91 = 4
        else:
            inel_mt_91 = 6
        rxn_91 = RXN.from_mf(mf=inel_mt_91, mt=91, mtname='inelastic(continuum)')
        return [rxn_min, rxn_max, rxn_91]

    else:

        return [rxn_min, rxn_max]
    

