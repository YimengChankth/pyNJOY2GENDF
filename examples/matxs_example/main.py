"""
Minimal example demonstrating macroscopic cross-section extraction using MatXS.

This example shows how to:
1. Generate GENDF files from ENDF data using NJOY
2. Define materials with MaterialSpec
3. Extract macroscopic cross-sections from GENDF files
4. Generate output CSV files for reactor physics calculations
"""

import pyNJOY2GENDF.njoy2gendf
from pyNJOY2GENDF.MatXS import MatXS, MaterialSpec
from pathlib import Path
import shutil
import sys


def find_njoy_executable():
    """Search for NJOY21 executable in common locations."""
    
    # List of potential paths to check
    potential_paths = [
        Path.home() / 'NJOY21' / 'bin' / 'njoy21',  # User home directory
        Path.home() / 'NJOY21' / 'build' / 'njoy21',  # User home directory
        Path.home() / 'bin' / 'njoy21',  # User bin
    ]
    
    # Check hardcoded and common paths first
    for path in potential_paths:
        if Path(path).is_file() and Path(path).stat().st_mode & 0o111: # check if file has execute permissions
            return str(path)
    
    # Check if njoy21 is in PATH
    njoy_in_path = shutil.which('njoy21')
    if njoy_in_path:
        return njoy_in_path
    
    return None


def generate_gendf_files(njoy_exec=None):
    """Generate GENDF files from ENDF data using NJOY."""
    
    print("=" * 60)
    print("Step 1: Generating GENDF Files")
    print("=" * 60)
    print("This will process ENDF data with NJOY to create GENDF files...")
    print("This may take several minutes depending on the number of nuclides.\n")
    
    # Find NJOY executable if not provided
    if njoy_exec is None:
        njoy_exec = find_njoy_executable()
        if njoy_exec is None:
            print("ERROR: NJOY21 executable not found!")
            print("\nPlease install NJOY21 or specify the path to the executable.")
            print("You can:")
            print("  1. Install NJOY21 and add it to your PATH")
            print("  2. Edit this script and set njoy_exec parameter in generate_gendf_files()")
            print("  3. Download from: https://github.com/njoy/NJOY21")
            sys.exit(1)
        else:
            print(f"Found NJOY executable: {njoy_exec}\n")
    else:
        # Verify provided path exists
        if not Path(njoy_exec).is_file():
            print(f"ERROR: NJOY executable not found at: {njoy_exec}")
            sys.exit(1)
        print(f"Using NJOY executable: {njoy_exec}\n")
    
    # Initialize the GENDF generator
    wr = pyNJOY2GENDF.njoy2gendf.njoy2gendf(
        savedir='gendf_library',                                         # Local directory to save files
        path2endf='../endfb8_example/ENDFB8_data/neutrons',          # Path to ENDF data
        njoy_exec=njoy_exec,                                             # NJOY executable
        nuclide_list=['Am241', 'B10'],  # Nuclides needed for examples
        filename_convention='ENDF-B'
    )
    
    # Define energy structure and conditions
    energy_bin_edges = [1e-5, 1, 100, 1000, 2e6]
    temperatures = [300, 600, 900]
    sig0s = [1e10, 1000, 1]
    
    # Generate NJOY input files and run NJOY
    wr.write_njoy_gendf_inputs(
        energybin_edges=energy_bin_edges,
        library_name='ENDF-B-VIII.0',
        temperatures=temperatures,
        sig0s=sig0s,
    )
    wr.run_njoy_inputs()
    
    print("\n✓ GENDF files generated successfully!")
    print(f"  Location: gendf_library/GENDF/\n")


def example_single_isotope_material():
    """Extract macro XS for a single-isotope material Am-241"""
    
    print("=" * 60)
    print("Example 1: Single Isotope Material")
    print("=" * 60)
    
    # Define a simple material with one isotope
    material = MaterialSpec(
        name='Am241_PURE',
        nuclides=['Am241.gendf'],
        number_densities=[4.8e-2],  # atoms/barn-cm
        temperature=300.0,  # Kelvin
        gendf_directory=Path('gendf_library/GENDF')
    )
    
    # Create processor and extract
    processor = MatXS([material], printmicro=False)
    results = processor.process_to_csv(output_directory='output_single')
    
    print("\n✓ Output written to: output_single/material_xs/Am241_PURE/")
    print("  Files: total.csv, fission.csv, scatter.csv, etc.")
    print(f"  Total XS (first 3 groups): {results['Am241_PURE']['macro']['total'][:3]}\n")


def example_multi_isotope_material():
    """Extract macro XS for a multi-isotope fuel material."""
    
    print("=" * 60)
    print("Example 2: Multi-Isotope Fuel")
    print("=" * 60)
    
    # Define fuel with Am241 and B10
    fuel = MaterialSpec(
        name='B10Am241_FUEL',
        nuclides=['Am241.gendf', 'B10.gendf'],
        number_densities=[
            4.5e-3,   # Am-241
            2.0e-2,   # B-10
        ],
        temperature=900.0,  # Fuel temperature
        gendf_directory=Path('gendf_library/GENDF')
    )
    
    processor = MatXS([fuel], printmicro=False)
    results = processor.process_to_csv(output_directory='output_multiple')
    
    print("\n✓ Output written to: output_multiple/material_xs/B10Am241_FUEL/")
    print("  - Includes sigma0 iteration for multi-isotope self-shielding")
    print("  - isolist.csv shows composition and temperature")
    print(f"  - Computed for {len(results['B10Am241_FUEL']['macro']['total'])} energy groups\n")


def example_multiple_materials():
    """Extract macro XS for multiple materials in one run."""
    
    print("=" * 60)
    print("Example 3: Multiple Materials")
    print("=" * 60)
    
    gendf_dir = Path('gendf_library/GENDF')
    
    # Define fuel
    fuel = MaterialSpec(
        name='FUEL',
        nuclides=['Am241.gendf', 'B10.gendf'],
        number_densities=[4.5e-3, 2.0e-2],
        temperature=900.0,
        gendf_directory=gendf_dir
    )
    
    # Define coolant
    coolant = MaterialSpec(
        name='COOLANT',
        nuclides=['B10.gendf'],
        number_densities=[2.2e-2],
        temperature=600.0,
        gendf_directory=gendf_dir
    )
    
    # Process all materials at once
    processor = MatXS([fuel, coolant], printmicro=False)
    results = processor.process_to_csv(output_directory='output_multi')
    
    print("\n✓ Output written to: output_multi/material_xs/")
    print(f"  Materials: {', '.join(results.keys())}")
    print("  - Each material in separate subdirectory")
    print("  - material_list file contains all material names\n")


def example_with_microscopic_output():
    """Extract both macro and micro XS."""
    
    print("=" * 60)
    print("Example 4: Macro + Microscopic XS")
    print("=" * 60)
    
    material = MaterialSpec(
        name='B10Am241_WITH_MICRO',
        nuclides=['Am241.gendf', 'B10.gendf'],
        number_densities=[4.5e-3, 2.0e-2],
        temperature=600.0,
        gendf_directory=Path('gendf_library/GENDF')
    )
    
    # Enable microscopic output with printmicro=True
    processor = MatXS([material], printmicro=True)
    results = processor.process_to_csv(output_directory='output_micro')
    
    print("\n✓ Output written to: output_micro/material_xs/B10Am241_WITH_MICRO/")
    print("  - macro/ directory: macroscopic cross-sections")
    if 'micro' in results['B10Am241_WITH_MICRO']:
        micro_keys = list(results['B10Am241_WITH_MICRO']['micro'].keys())
        print(f"  - Microscopic data for: {', '.join(micro_keys)}")
    print("  - Results also available in returned dictionary\n")


def example_from_monte_input():
    """Use legacy MONTE input file format (if available)."""
    
    print("=" * 60)
    print("Example 5: From MONTE Input File")
    print("=" * 60)
    
    # This example assumes you have a MONTE input file
    # Format:
    #   XSP "path/to/GENDF"
    #   SIGNAL TEMP CONSTANT 900.0
    #   MAT FUEL U235.gendf 4.5e-3 TEMP
    #   MAT FUEL U238.gendf 2.0e-2 TEMP
    
    input_file = Path('monte_input.inp')
    
    if input_file.exists():
        processor = MatXS.from_monte_input(input_file, printmicro=False)
        results = processor.process_to_csv(output_directory='output_monte')
        print(f"\n✓ Successfully processed MONTE input file")
        print(f"  Materials processed: {', '.join(results.keys())}\n")
    else:
        print(f"\n⚠ MONTE input file not found: {input_file}")
        print("  This example requires a MONTE-format input file.")
        print("  See MONTE input documentation for format details.\n")


if __name__ == '__main__':
    print("\n" + "=" * 60)
    print("MatXS: Complete GENDF Generation and XS Extraction")
    print("=" * 60 + "\n")
    
    # Step 1: Generate GENDF files (comment this out if files already exist)
    # If you need to specify a custom NJOY path, use:
    # generate_gendf_files(njoy_exec='/path/to/your/njoy21')
    generate_gendf_files()
    
    # Step 2: Run MatXS examples (uncomment the ones you want to run)
    print("\n" + "=" * 60)
    print("Step 2: Extracting Macroscopic Cross-Sections")
    print("=" * 60 + "\n")
    
    # Uncomment the examples you want to run:
    # example_single_isotope_material()
    example_multi_isotope_material()
    # example_multiple_materials()
    # example_with_microscopic_output()
    # example_from_monte_input()
    
    print("\n" + "=" * 60)
    print("Instructions:")
    print("=" * 60)
    print("1. Uncomment the desired example function(s) above")
    print("2. Run: python3 main.py")
    print("3. Check the output_* directories for results")
    print("\nNote: If GENDF files already exist, comment out the")
    print("      generate_gendf_files() call to save time.")
    print("\nTo specify a custom NJOY path:")
    print("      generate_gendf_files(njoy_exec='/path/to/njoy21')")
    print("=" * 60 + "\n")
