from pyNJOY2GENDF.ioreadwrite import read_endf_mt_mf, print_tuples
import numpy as np 

class WeightingFlux:
    '''Weighting flux for GROUPR card 8b see page 238
        • An ENDF/B “TAB1” record consists of three distinct parts:
        two double values and four integer values of which only the last two integers
        (the number of interpolation ranges NR and the number of E,C(E) pairs
        NP) are needed to read the remainder of the record

        • the interpolation scheme data which is a sequence of NR pairs of NBT (index
        of the E,C(E) pair corresponding to the end of interpolation range) and
        INT (the interpolation type)

        • the tabulated data which is a sequence of NP E,C(E) pairs

        0. 0. 0 0 NR NP
        NBT(1) INT(1) ... NBT(NR) INT(NR)
        E(1) C(1) ... E(NP) C(NP) / card8b

    '''

    def __init__(self, interpolate_range:list[tuple[float, float]] = [], flux_vals:list[tuple[float, float]] = []):
        '''
        Docstring for __init__
        
        :param self: Description
        :param interpolate_ints: list of (NBT, INT) values. Note this should be the exact values seen in the manual. NBT should be cumulative index not the # of energy groups. See also compress_vector 
        :param flux_vals: Description : list of (E, C(E)) values. E(1) should be the lowest energy!
        '''

        # hint: int values can take on the following
        # 1 : histogram
        # 2 : linear-linear
        # 3 : linear-log
        # 4 : log-linear  
        # 5 : log-log (if there is doubt, use 5)

        self.interpolate_range = interpolate_range
        self.flux_vals = flux_vals

        pass
    
    @staticmethod   
    def compress_vector(v: list[int]) -> list[tuple[int, int]]:
        """
        Compress a vector into (end_index, value) pairs using 1-based indexing.

        Example
        -------
        [2, 2, 2, 1, 1, 1, 1] -> [(3, 2), (7, 1)]
        """

        if not v:
            return []

        result: list[tuple[int, int]] = []

        current_value = v[0]
        count = 1

        for x in v[1:]:
            if x == current_value:
                count += 1
            else:
                result.append((count + (result[-1][0] if result else 0), current_value))
                current_value = x
                count = 1

        # add last run
        result.append((count + (result[-1][0] if result else 0), current_value))

        return result


    def write_to_block(self):
        NR = len(self.interpolate_range) # number of interpolation ranges 
        NP = len(self.flux_vals) # Number of E,C(E) pairs
        lines = f'   0.0 0.0 0 0 {NR} {NP}\n'
        # interpolation specification,     
        lines += print_tuples(self.interpolate_range, break_every=3, rjust_size=12)
        if (lines.endswith('\n') == False):
            lines = lines + '\n'  # ensure line break
        flux_vals_scientific = [(f'{e:1.5e}', f'{c:1.5e}') for e, c in self.flux_vals]
        lines += print_tuples(flux_vals_scientific, break_every=3, rjust_size=12)
        lines = lines + ' / card8b\n'    
        # ensure line break
        return lines




class GROUPRCard8:
    def __init__(self):
        self.fehi_dict = None
        self.sigpot_dict = None
        self.nflmax = None
        self.weightingflux = None
        self.iwt = None
        self.alpha2 = 0.0
        self.sam = 0.0
        self.beta = 0.0
        self.alpha3 = 0.0
        self.gamma = 0.0
        self.jsigz = 0
        self.ninwt = 0

        pass

    @property
    def iwt(self):
        return self._iwt
    @iwt.setter
    def iwt(self, value):
        self._iwt = value

    @property
    def fehi(self):
        '''
        Only useful for iwt = -n (card 8a)
        
        :param self: Description
        '''
        return self._fehi
    
    @fehi.setter
    def fehi(self, value):
        self._fehi = value

    @property
    def sigpot(self):
        '''
        Only useful for iwt = -1 (card 8a)
        
        :param self: Description
        '''
        return self._sigpot
    @sigpot.setter
    def sigpot(self, value):
        self._sigpot = value

    @property
    def weightingflux(self):
        return self._weightingflux
    @weightingflux.setter
    def weightingflux(self, value): 
        assert isinstance(value, (type(None), WeightingFlux)), 'weightingflux must be WeightingFlux object or None'
        self._weightingflux = value


    @classmethod
    def from_endf(cls, endfdata:str, nflmax=20000) -> tuple['GROUPRCard8', int]:
        '''

        IMPORTANT: This is very developmen

        Create a GROUPRCard8 object with default values from ENDF data file

        Looking at https://www-nds.iaea.org/wimsd/processing.htm, it seems like only fission products and actinides use Card 8 weighting flux. 


        '''
        raise NotImplementedError('GROUPRCard8.from_endf method is under development and probably will not be used in future. Use manual input for now.')
        print(f'Reading ENDF data to determine Card 8. Path:{endfdata}', flush=True)

        mt151 = read_endf_mt_mf(endfdata, MT=151, MF=2)  # check if resonance data exists. If it does not exist, then use iwt = 1 (Bodorenko flux), else -1
        # see pg 447 of NJOY manual for description on how to find upper range of resolved resonance and scattering potential
        # following WIMS suggested input files at \url{https://www-nds.iaea.org/wimsd/processing.htm}, assume that for nuclei without resonance scattering, use narrow resonance approximation Bodorenko flux (iwt=1)

        # bodorenko flux calc \sigma(\sigma_{b})
        # \phi(E) = C(E)/(\sigma_{t}(E) + \sigma_{0}) where \sigma_{0} is the background cross-section

        # Flux weighting calculation  (see pages 197-200) for E < fehi
        # homogeneous flux problem is given by eqn 272 on page 198
        # [\sigma_{0} + \sigma_{t}(E)]\phi(E) = C(E)\sigma_{0} + \int_{E}^{E/\alpha} \frac{\sigma_{s}(E')}{(1-\alpha)E'}\phi(E')dE'
        # This requires the following parameters:
        # fehi : break between computed flux and bondarenko flux (must be in the resolved resonance region)
        # sigpot : potential scattering cross-section in barns. "The actual value for sigpot is not very critical - a number near 10 barns is typical for fissionable materials."
        # nflmax : maximum number of flux iterations (suggested value 5000 in WIMS)

        assert mt151 is not None, 'ENDF data must contain MT151 MF2 to use GROUPRCard8.from_endf method'

        fehi = mt151[2, 1]  # second value shows the upper range of the resolved resonance energy fehi
        # for many nuclides, the upper range can be very large. e.g. Pu-239 it is 2500 eV. But in WIMS file fehi is set to 622.2. So where does this 622.2 eV come from??

        # case 1: fehi < 600 (typical of heavy nuclei)
        # see section 8.5


        fehi_threshold = 1000  # eV

        if fehi < fehi_threshold:
            print(f'\nMT=151,MF=2 has fehi < {fehi_threshold}, treatment using flux calculator method with tabulated user supplied flux (iwt=-1)', flush=True)
            r = mt151[3, 1] # scattering potential radius. 
            sigpot = 4*np.pi*(r**2) # see pg 447 NJOY manual for formula. "The actual value for sigpot is not very critical - a number near 10 barns is typical for fissionable materials."
            weightingflux = None  # Create a None object and let user fill in later
            gprc8 = GROUPRCard8()
            gprc8.fehi = fehi
            gprc8.sigpot = sigpot
            gprc8.nflmax = nflmax
            gprc8.weightingflux = weightingflux
            gprc8.iwt = -1
            print(f'Params: \n\tfehi: {fehi}, \n\tsigpot: {sigpot:1.6e}, \n\tnflmax: {nflmax}', flush=True)
            return gprc8  # iwt = -1 indicates the use of flux calculator method in NJOY
        else:
            print(f'MT=151,MF=2 has fehi:{fehi} >= {fehi_threshold}, treatment using bodorenko flux method with tabulated user supplied flux (iwt=1)', flush=True)
            # case 2: bodorenko flux calc (typical of light nuclei). See
            # section 8.4. Then Card8a is not needed and only card  weighting flux is needed 
            gprc8 = GROUPRCard8()
            gprc8.iwt = 1
            return gprc8  # iwt = 1 indicates the use of bodorenko flux calculator method in NJOY


    def write_to_block(self):
        # check if card 8a is needed 
        # Cards 1-7 of GROUPR is handled by templates.GROUPR
        tmp = ''

        if self.iwt <= -1:  # if iwt is negative, then card 8a is needed. If iwt is positive, then card 8a is not needed. If iwt = 0, then neither card 8a nor 8b is needed:
            tmp += f'   {self.fehi} {self.sigpot} {self.nflmax} {self.ninwt} {self.jsigz} {self.alpha2} {self.sam} {self.beta} {self.alpha3} {self.gamma}  / card8a fehi, sigpot, nflmax, ninwt, jsigz, alpha2, sam, beta, alpha3, gamma \n'
        if self.iwt in [1, -1]:
            tmp += self.weightingflux.write_to_block() # this is card 8b, weighting flux
        return tmp

