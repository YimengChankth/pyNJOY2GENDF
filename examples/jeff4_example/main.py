import pyNJOY2GENDF.njoy2gendf

if __name__ == "__main__":

    wr = pyNJOY2GENDF.njoy2gendf.njoy2gendf(savedir='jeff4_gendf_library', # path to save the files to 
                    path2endf='JEFF4_data',                                # folder of raw data 
                    njoy_exec='/data/user/chan_y/NJOY21/bin/njoy21',       # change the njoy21 executable here! You can install from github 
                    nuclide_list=['Am241','B10'],                          # list the nuclides you want to process
                    filename_convention='JEFF')
    
    energy_bin_edges = [1e-5, 1, 100, 1000, 2e6]
    temperatures = [300, 600, 900]
    sig0s = [1.0]

    wr.write_njoy_gendf_inputs(
        energybin_edges=energy_bin_edges,
        library_name='JEFF4.0', # this is cosmetic i.e. not relevant to calcualtions 
        temperatures=temperatures,
        sig0s=sig0s,
    )
    wr.run_njoy_inputs()