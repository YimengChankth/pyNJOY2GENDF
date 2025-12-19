from abc import ABC, abstractmethod
import numpy as np

# abstract class for RECONR
class RECONR(ABC):
    @abstractmethod
    def write_block(self):
        pass

# default for RECONR (based on KM)
class default_RECONR(RECONR):
    @staticmethod
    def write_block(nuclide='Am241', mat=9543, library='ENDF-B/VIII'):
        tmp =f'''-- 1 -- -- 2 -- -- 3 -- -- 4 -- -- 5 -- -- 6 -- -- 7 -- -- 8 -- -- 9 -- -- 0 --
RECONR
   20      21 / nendf   npend
   '{nuclide} {library}' / tlabel
   {str(mat)} / mat
   0.01 / err
   0 / end\n
'''
        
        return tmp
        
# abstract class for BROADR
class BROADR(ABC):
    @abstractmethod
    def write_block(self, temperatures):
        pass

# default for RECONR (based on KM)
class default_BROADR(BROADR):
    @staticmethod
    def write_block(mat, temperatures):
        ntemp = len(temperatures)

        temp_str = '    '.join([f'{i:.1f}' for i in temperatures])

        tmp = f'''-- 1 -- -- 2 -- -- 3 -- -- 4 -- -- 5 -- -- 6 -- -- 7 -- -- 8 -- -- 9 -- -- 0 --
BROADR
   20      21      22 / nendf nin nout
   {str(mat)}    {ntemp} / mat1 ntemp1
   0.01 / errthn
   {temp_str} / temp2
   0 / end\n
'''
        return tmp
        
# abstract class for UNRESR     
class UNRESR(ABC):
    @abstractmethod
    def write_block(self, temperatures, sig0s):
        pass

class default_UNRESR(UNRESR):
    @staticmethod
    def write_block(mat, temperatures, sig0s):
        ntemp = len(temperatures)
        nsig0s = len(sig0s)

        temp_str = '    '.join([f'{i:0.1f}' for i in temperatures])
        sig0s_str = '    '.join([f'{i:.3e}' for i in sig0s])

        tmp = f'''-- 1 -- -- 2 -- -- 3 -- -- 4 -- -- 5 -- -- 6 -- -- 7 -- -- 8 -- -- 9 -- -- 0 --
UNRESR
   20      22      23 / nendf nin nout
   {str(mat)}    {ntemp}  {nsig0s} / matd ntemp nsigz
   {temp_str} / temp
   {sig0s_str} / sigz
   0 / end\n
'''

        return tmp


class HEATR(ABC):
    @abstractmethod
    def write_block(self):
        pass

# I changed from KM's definition, since we dont actually need the fission kerma, we can just use the total kerma (mt=301)
class default_HEATR(HEATR):
    @staticmethod
    def write_block(mat):
        tmp = f'''-- 1 -- -- 2 -- -- 3 -- -- 4 -- -- 5 -- -- 6 -- -- 7 -- -- 8 -- -- 9 -- -- 0 --
HEATR
   20      23      24 / nendf nin nout
   {str(mat)}    0 / matd npk\n
'''
        return tmp 
    
# abstract base class for GROUPR which is the most complicated 
class GROUPR(ABC):
    @abstractmethod
    def write_block(self, mat:str, library:str, temperatures, sig0s, group_edges):
        pass




class default_GROUPR(GROUPR):
    @staticmethod
    def write_block(nuclide:str, mat:int, library:str, temperatures:list[float], sig0s:list[float], energybin_edges:list[float], rxn_base:list['RXN'], rxn_temp:list['RXN']):
        '''

        Example of output

        -- 1 -- -- 2 -- -- 3 -- -- 4 -- -- 5 -- -- 6 -- -- 7 -- -- 8 -- -- 9 -- -- 0 --
        GROUPR
           20      24      0       25 / nendf npend ngout1 ngout2
           9543    1       0       8       8       10      7 / matb ign igg iwt lord ntemp nsigz
           'Am241 JEFF3.3' / tlabel
           300.    600.    900.    1200.   1500.   1800.   2100.   2400.   2700.   3000. / temp
           1.e10   1.e5    1.e4    1000.   100.    10.     1. / sigz
           25 / ngn
           2.15443E-1      4.64159E-1      1.00000E+0      2.15443E+0      4.64159E+0      1.00000E+1
           2.15443E+1      4.64159E+1      1.00000E+2      2.15443E+2      4.64159E+2      1.00000E+3
           2.15443E+3      4.64159E+3      1.00000E+4      2.15443E+4      4.64159E+4      1.00000E+5
           2.00000E+5      4.00000E+5      8.00000E+5      1.40000E+6      2.50000E+6      4.00000E+6
           6.50000E+6      10.5000E+6 / egg
           3       1       total / mfd mtd mtname
           3       102     ngamma /
           3       18      fission /
           8       16      n2n / 
           6       2       elastic /
           6       51      inelastic /
           6       -90     inelastic /
           8  91  inelastic  /
           3       259     inv /
           3       452     nubar /
           5       18      chi /
           3       301     heat /


        Parameters
        ----------
        nuclide : string 
            isotope. name cosmetic 

        mat : int
            mat id of isotope

        library : str
            cosmetic name

        temperatures : list of floats
            temperature to evalaute 

        sig0s : list of floats
            sig0s to evalaute 
        
        group_edges : list of floats
            energy group bin edges in eV
        
        rxn_base : list of ints 
            reaction MT numbers for first temperature 

        rxn_temp : list of ints 
            reaction MT numbers for other temperatures. Normally you do temperature dependence for a small subset of MTs 
        '''
        ntemp = len(temperatures)
        nsig0s = len(sig0s)
        ngn = len(energybin_edges) - 1
        temp_str = '    '.join([f'{i:0.1f}' for i in temperatures])
        sig0s_str = '    '.join([f'{i:.3e}' for i in sig0s])

        header=f'''-- 1 -- -- 2 -- -- 3 -- -- 4 -- -- 5 -- -- 6 -- -- 7 -- -- 8 -- -- 9 -- -- 0 --
GROUPR
   20      24      0       25 / nendf npend ngout1 ngout2
   {str(mat)}    1       0       8       8    {ntemp}       {nsig0s}     / matb ign igg iwt lord ntemp nsigz
   '{nuclide} {library}' / tlabel
   {temp_str} / temp
   {sig0s_str} / sigz
   {ngn} / ngn
'''
        group_edges_str = '   '
    
        for cge, ge in enumerate(energybin_edges):
            group_edges_str += f'{ge:1.5e}'
            if (cge % 6) != 5:
                group_edges_str += ' '*6
            else:
                group_edges_str += '\n   '
        group_edges_str += '/ egg\n'

        rxn_string = ''
        # reactions for base temperature 
        for rxn in rxn_base:
            rxn_string += '   ' + rxn.write_line() + '\n'
        rxn_string +='   0 /\n'

        # reactions for nonbase temperature
        for it in range(ntemp - 1):
            for rxn in rxn_temp:
                rxn_string += '   ' + rxn.write_line() + '\n'
            
            rxn_string +='   0 /\n'

        return header + group_edges_str + rxn_string + '   0 / end\n'


class RXN:
    def __init__(self, mfd, mt, mtname):
        self.mfd = mfd
        self.mt = mt
        self.mtname = mtname

    @classmethod
    def from_mf(cls, mf, mt, mtname) -> 'RXN':
        r'''
        # Elastic and inelastic scattering explanation for MT=4/5 -> MT=6
        MT = 4 is the energy independent scattering distribution of secondary particle
        MT = 5 is the scattering angle independent distribution of energy

        With assumption of independence, the joint distribution is calculated in NJOY using MFD=6 (see Eqn 328, page 218 in NJOY manual). This requires that MT=[4,5] exists in ENDF file for MT of interest. 
        $$
        f(E\rightarrow E', \mu) = F_{4}(E, \mu)F_{5}(E\rightarrow E')
        $$

        For some isotopes, the joint probability is given directly in MT=6. 
        MT = 6 is the joint scattering distribution and energy distribution i.e. 
        $$
        f(E\rightarrow E', \mu) = F_{4}(E, \mu)F_{5}(E\rightarrow E') = F_{6}(E\rightarrow E', \mu)
        $$
        NJOY processes this using MFD=8, which requires that MT=6 exists in ENDF file for MT number. Note that elastic scattering will always have MT=6 (energy distribution out is ignored).

        '''

        mf2mfd = {
            3:3,
            4:6,
            5:6,
            6:8,
        }

        mfd = mf2mfd[mf]
        return RXN(mfd, mt, mtname)
    
    def write_line(self):
        return f'{self.mfd}'.ljust(8) + f'{self.mt}'.ljust(8) + self.mtname + ' /'


def mfd_describ():
    '''
    A handy guide to MFD numberes meaning from the NJOY handbook. See page 244 in manual. 
    '''

    mfd_des = {
        6:'joint energy out, angle distribution of secondary particle',
        259:'ave. inverse neutron velocity for group in s/m',
        452:'average total fission yield computed from MF=1 and MF=3',
        301:'heat production'
    }

    pass

