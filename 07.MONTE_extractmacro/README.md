To compile on Merlin7:
module load PrgEnv-gnu
make 

1) Copy of MONTE but only in the input file reading and export cross-sections to csv file for processing
2) Python script to read the relavant csv files and to create openmc.mgxs libraries 

## scattering problem

In MONTE_extractmacro.f90, I print out the scattering matrices: scatter_0 is the 0-order expansion, scatter_1 is the 1st order expansion etc. The scattering matrices are the same.

## ym notes
OpenMC multi-group notes: 

Ref: https://docs.openmc.org/en/v0.10.0/examples/mg-mode-part-i.html 


### Create a 7-group structure with arbitrary boundaries (the specific boundaries are unimportant)
groups = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, 8))

uo2_xsdata = openmc.XSdata('uo2', groups)
uo2_xsdata.order = 0

### When setting the data let the object know you are setting the data for a temperature of 294K.
uo2_xsdata.set_total([1.77949E-1, 3.29805E-1, 4.80388E-1, 5.54367E-1,
                      3.11801E-1, 3.95168E-1, 5.64406E-1], temperature=294.)

uo2_xsdata.set_absorption([8.0248E-03, 3.7174E-3, 2.6769E-2, 9.6236E-2,
                           3.0020E-02, 1.1126E-1, 2.8278E-1], temperature=294.)
uo2_xsdata.set_fission([7.21206E-3, 8.19301E-4, 6.45320E-3, 1.85648E-2,
                        1.78084E-2, 8.30348E-2, 2.16004E-1], temperature=294.)

uo2_xsdata.set_nu_fission([2.005998E-2, 2.027303E-3, 1.570599E-2, 4.518301E-2,
                           4.334208E-2, 2.020901E-1, 5.257105E-1], temperature=294.)

uo2_xsdata.set_chi([5.87910E-1, 4.11760E-1, 3.39060E-4, 1.17610E-7,
                    0.00000E-0, 0.00000E-0, 0.00000E-0], temperature=294.)

the shape of OpenMC's scattering matrix entry is instead [Incoming groups, Outgoing Groups, Scattering Order]
set_scatter_matrix






# Initialize the library
mg_cross_sections_file = openmc.MGXSLibrary(groups)

# Add the UO2 data to it
mg_cross_sections_file.add_xsdata(uo2_xsdata)

# And write to disk
mg_cross_sections_file.export_to_hdf5('mgxs.h5')
