import openmc

def run_openmc():
    # MATERIALS 
    materials = {}
    for xs in ['INNER','SHIELD']:
        materials[xs] = openmc.Material(name=xs)
        materials[xs].set_density('macro', 1.)
        materials[xs].add_macroscopic(xs)

    # Instantiate a Materials collection, register all Materials, and export to XML
    materials_file = openmc.Materials(materials.values())

    # Set the location of the cross sections file to our pre-written set
    materials_file.cross_sections = '../INP/material_xs/mgxs.h5'