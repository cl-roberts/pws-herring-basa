# TODO

Note: all file paths reported in this TODO are relative to `pws-herring-basa/sensitivity`.
Don't get confused with the base model code in `pws-herring-basa`!

## List of design tasks

1. [ ] Fill in user stories
2. [ ] Maybe some more implicit use cases?
3. [ ] Fill in components
4. [x] Create organizational design chart for `design.md` 

## List of package components to implement

Below is a list of `python` modules that we will need to write and save to 
`src/sensitivity`. Each module should create a function(s) of the same/similar name that 
has a docstring(s) and associated tests saved to `tests` with the naming convention 
`test_module.py`. The tests can 
be super basic. At the minimum, they should simply exist and test to see if a file 
was created (exceptions noted below). The important thing is to have the structure 
and workflow in place so that we have a clear roadmap for extending the codebase 
if we had more time to work on it.

1. [x] `run_basa.py` executes the model and returns clean output data. 
   - [ ] Test
2. [ ] `modify_mortality.py` changes assumed rate of natural mortality in model
   - [ ] Test
3. [ ] `plot_sensitivity.py` plots sensitivity model biomass time series alongside 
       base model biomass
   - [ ] Test
4. [ ] `read_biomass.py` wrappers to read biomass estimates from `.csv`
   - [ ] Test
5. [ ] `sensitivity_comparison.py` calculates percent error (or any other relevant 
       metric) between base and sensitivity biomass and reports summary statistics
   - [ ] Test
6. [ ] `plot_sensitivity_comparison.py` plots percent error between base and sensitivity
       biomass by year
   - [ ] Test


## List of shiny app components to implement

### User interface

1. [x] Natural mortality slider bar input
2. [ ] Natural mortality number box input
3. [ ] Other buttons/sliders/etc. for other parameters as place holders
4. [x] Boxes to report error metric summary statistics
5. [ ] Box to show sensitivity/base biomass time series plot
6. [x] Box to show error time series plot
7. [ ] Tabs to create multiple sensitivity runs
8. [ ] Action button to initiate sensitivity analysis
9. [ ] Download button to obtain sensitivity results

