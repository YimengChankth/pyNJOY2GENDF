# MatXS Example

This example demonstrates how to extract macroscopic cross-sections from GENDF files using the `MatXS` class.

## API Overview

The `MatXS` class provides two methods for processing materials:

1. **`process()`** - Returns cross-section data as Python dictionaries (no file I/O)
2. **`process_to_csv(output_directory)`** - Computes data and writes CSV files to disk

This separation allows you to:
- Use cross-section data directly in Python without writing files
- Integrate with other tools and workflows programmatically
- Still generate traditional CSV output when needed

## Prerequisites

Before running this example, you must first generate GENDF files using one of the other examples:

```bash
cd ../endfb8_example
python3 main.py
# Wait for NJOY to complete...
```

Or for JEFF-4.0:

```bash
cd ../jeff4_example
python3 main.py
```

## Running the Examples

Edit `main.py` and uncomment the example function(s) you want to run:

```python
# Uncomment any of these:
example_single_isotope_material()
example_multi_isotope_material()
example_multiple_materials()
example_with_microscopic_output()
example_from_monte_input()
```

Then execute:

```bash
python3 main.py
```

## Examples Included

### 1. Single Isotope Material
Extract macroscopic XS for a pure material (e.g., pure U-235).

### 2. Multi-Isotope Fuel
Extract macro XS for UO₂ fuel with uranium enrichment. Demonstrates sigma0 iteration for self-shielding in multi-isotope materials.

### 3. Multiple Materials
Process fuel, coolant, and cladding materials in a single run. Output organized by material name.

### 4. Macro + Microscopic XS
Generate both macroscopic (homogenized) and microscopic (per-isotope) cross-sections.

### 5. From MONTE Input File
Parse legacy MONTE input format (requires input file in MONTE format).

## Output Structure

### Using `process()` - Dictionary Output

The `process()` method returns a dictionary with this structure:

```python
{
    '<material_name>': {
        'name': str,              # Material name
        'temperature': float,      # Temperature in Kelvin
        'nuclides': list[str],    # List of nuclide filenames
        'number_densities': ndarray,  # Atom densities (atoms/barn-cm)
        'sig0': ndarray,          # Background cross-sections (niso, ng)
        'macro': {
            'total': ndarray,         # Total XS (ng,)
            'fission': ndarray,       # Fission XS (ng,)
            'nu_fission': ndarray,    # Nu*fission XS (ng,)
            'chi': ndarray,           # Fission spectrum (ng,)
            'kerma': ndarray,         # KERMA factors (ng,)
            'absorption': ndarray,    # Absorption XS (ng,)
            'scatter': ndarray        # Scattering matrix (nl, ng, ng)
        },
        'micro': {  # Only if printmicro=True
            '<nuclide_name>': {
                'total': ndarray,
                'fission': ndarray,
                'kerma': ndarray,
                'nu': ndarray,
                'chi': ndarray,
                'flux': ndarray  # Optional
            }, ...
        }
    }, ...
}
```

### Using `process_to_csv()` - File Output

Each example creates an output directory with structure:

```
output_<example>/
└── material_xs/
    ├── material_list
    └── <material_name>/
        ├── isolist.csv
        ├── sig0.csv
        └── macro/
            ├── total.csv
            ├── fission.csv
            ├── nu_fission.csv
            ├── chi.csv
            ├── absorption.csv
            ├── kerma.csv
            └── scatter.csv
```

## Usage Examples

### Working with Dictionaries (In-Memory)

```python
from pyNJOY2GENDF.MatXS import MatXS, MaterialSpec
from pathlib import Path

# Define material
material = MaterialSpec(
    name='UO2_FUEL',
    nuclides=['U235.gendf', 'U238.gendf', 'O16.gendf'],
    number_densities=[4.5e-3, 2.0e-2, 4.2e-2],
    temperature=900.0,
    gendf_directory=Path('gendf_library/GENDF')
)

# Process and get results as dictionaries
processor = MatXS([material], printmicro=False)
results = processor.process()

# Access data directly
fuel_data = results['UO2_FUEL']
total_xs = fuel_data['macro']['total']
chi = fuel_data['macro']['chi']
scatter_matrix = fuel_data['macro']['scatter']
```

### Writing to CSV Files

```python
# Process and write CSV files
processor = MatXS([material], printmicro=True)
results = processor.process_to_csv(output_directory='output_fuel')

# Results are also returned as dictionaries
# So you can use them immediately without reading files
print(f"Total XS: {results['UO2_FUEL']['macro']['total']}")
```

## Customization

Modify the `MaterialSpec` parameters to suit your needs:
- `nuclides`: GENDF filenames
- `number_densities`: Atom densities in atoms/barn-cm
- `temperature`: Material temperature in Kelvin
- `gendf_directory`: Path to GENDF files

## Notes

- **New in this version**: `process()` returns data as dictionaries, `process_to_csv()` writes CSV files
- Sigma0 iteration is performed automatically for multi-isotope materials
- Temperature interpolation uses linear interpolation
- Sigma0 interpolation uses log10 interpolation
- Scattering matrices are written in sparse format to save disk space
- Use `printmicro=True` to get per-nuclide contributions in both dictionary and CSV outputs
