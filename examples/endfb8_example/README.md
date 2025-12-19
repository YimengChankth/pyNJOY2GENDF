To run this example to generate GENDF files for Am-241 and B-10 with 5 energy group structures with bins [1e-5, 1, 100, 1000, 2e6] ev at temperatures [300, 600]

!! Remember to change the NJOY21 executable path!! 

```
python3 main.py 
```

Expected outputs
- GENDF files located in folder 'endfb8_gendf_library'

    + Collection of desired GENDF files:
        -- endfb8_gendf_library/GENDF/Am241.gendf
        -- endfb8_gendf_library/GENDF/B10.gendf

    + Collection of NJOY inputs (check this to see if NJOY fails):
        -- endfb8_gendf_library/NJOY_INPUTS/Am241.inp [Notice that fission is automatically detected, so MT18 is included]
        -- endfb8_gendf_library/NJOY_INPUTS/B10.inp

    +  Collection of NJOY outputs:
        -- endfb8_gendf_library/NJOY_OUTPUTS/Am241.out
        -- endfb8_gendf_library/NJOY_OUTPUTS/B10.out

    + Summary of outputs (quick check if NJOY21 has reported any errors.)
        -- endfb8_gendf_library/NJOY_OUTPUTS/njoy_run_analysis.tex


