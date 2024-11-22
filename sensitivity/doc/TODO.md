# TODO

Note: all file paths reported in this TODO are relative to `pws-herring-basa/sensitivity`.
Don't get confused with the base model code in `pws-herring-basa`!

## List of design tasks

1. [ ] Fill in user stories
2. [ ] Maybe some more implicit use cases?
3. [ ] Fill in components

## List of package components to implement

Below is a list of `python` modules that we will need to write and save to 
`src/sensitivity`. Each module should create a function(s) of the same/similar name that 
has a docstring(s) and associated tests saved to `tests` with the naming convention 
`test_module.py`. The tests can 
be super basic. At the minimum, they should simply exist and test to see if a file 
was created (exceptions noted below). The important thing is to have the structure 
and workflow in place so that we have a clear roadmap for extending the codebase 
if we had more time to work on it.

1. [ ] `run_basa.py` executes the model and returns clean output data. 
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


### Module 1: `run_basa.py`

#### Executes the model and returns clean output data

This function will just be a wrapper that takes no arguments and executes that 
following shell commands:

```
Rscript run_basa.r
Rscript plotting/plot_management_outputs.r
```

This [stackoverflow thread](https://stackoverflow.com/questions/19894365/running-r-script-from-python)
may be relevant. 

I think that `run_basa()` should create a new directory (e.g. `data_outputs/sensitivity_1`) 
when executed and copy the sensitivity outputs to that directory. That way we can
do multiple sensitivity runs and cache the results of earlier runs. 

**Important note:** for some unusual values of assumed natural mortality, the model 
will fail to converge. This module should be able to identify when this is the case
and clearly communicate it to the user. It will probably need to run the commands
above in a `try` block. Perhaps this module should return a `0` if the model 
converges and a `1` otherwise? The test for this module should include a edge 
test that ensures the model run fails gracefully if it cannot fit to the mortality 
values it is given.  


### Module 2: `modify_mortality.py`

#### changes assumed rate of natural mortality in model

This module is makes the key changes to the model assumptions which provide 
the basis for our sensitivity analysis. Despite its importance, its scope will be 
very narrow and hopefully easy to implement. The natural mortality assumption
gets fed into the model via the file `model/PWS_ASA(par).ctl`. Note that a `.ctl`
file is basically a `.txt` file that ADMB uses to control model execution (we call
it a "control" file). Lines 12-19 in `model/PWS_ASA(par).ctl` look like this:

```
## —————————————————————————————————————————————————————————————————————————— ##
##  init   lower   upper    est  prior                fun                     ##
## value   bound   bound    phz   type     p1    p2  type   # PARAMETER       ##
## —————————————————————————————————————————————————————————————————————————— ##
    0.08    0.00    5.00      3      0      0     0     0   # 1:  VHSV_age3_4_mort_93
    0.22    0.00    5.00      3      0      0     0     0   # 2:  ICH_age5_8_mort_93
    0.25    0.05    1.50     -2      0      0     0     0   # 3:  Z_0_8
    0.93    0.30    1.60      2      0      0     0     0   # 4:  Z_9
```

The natural mortality parameter we want to modify is called `Z_0_8` in the model.
Note that I'm using $M$ to denote this parameter but the model calls it `Z_0_8`.
The other parameters shown above are other mortality parameters that are estimated
and we don't need to touch them. In the `est phz` column in the table above, a $-$
sign denotes a "fixed" parameter, so this control file is telling the model to 
fix $M$ at 0.25. 

So, we want `modify_mortality.py` to take `model/PWS_ASA(par).ctl` and a new value
for `Z_0_8` as inputs, correctly identify the $M =  0.25$ in the table above and 
change it to the new value, then write over `model/PWS_ASA(par).ctl` with the alteration.

#### A note on natural mortality:

In fisheries population models, natural mortality refers to the instantaneous 
background rate of death in the stock and is generally denoted $M$. We usually 
estimate or define $M$ for specific years ($M_t$) and can re-express it as an 
annual survival fraction, $S_t$, which gives the proportion of fish
that survived year $t$:

$$ S_t = e^{-M_t} $$

In BASA, $M_t$ is constant over time so the subscript is not relevant to us.
Since $S_t$ is a proportion, the $M_t$ is only valid for $M_t \in [0, \infty)$. 
So we should probably raise a `valueError` if `modify_mortality()` receives a 
negative number, but in practice the range of values for $M$ that users will be 
able to choose in the shiny app will be much narrower than $[0, \infty)$.

### Module 3. `read_biomass.py`

#### wrappers to read biomass estimates from `.csv`

This module should create two functions:

`read_biomass_base()`: takes no input, reads `data_outputs/outputs-for-management_base.csv`
and saves it to a pandas DataFrame.

`read_biomass_sensitivity()`: takes a directory input created by `run_basa()` that 
points to a particular sensitivity run (e.g. `data_outputs/sensitivity_1`), 
reads `outputs-for-management.csv` and saves to a pandas DataFrame. 

### Module 4. `plot_sensitivity.py`

#### plots sensitivity model biomass time series alongside base model biomass

This module will create a function that takes a directory input created by 
`run_basa()` that points to a particular sensitivity run (e.g. 
`data_outputs/sensitivity_1`) and use it to call `read_biomass_sensitivity()`.
Then it will call `read_biomass_base()` and plot both time series with
different colors and a legend. `plot_sensitivity()` will either return the plot or 
save it to a file (not sure which would be better).

### Module 5. `sensitivity_comparison.py`

#### calculates percent error (or any other relevant metric) between base and sensitivity biomass and reports summary statistics

This module will create a function which calls `read_biomass_sensitivity()` and
`read_biomass_base()` and calculate a percent error (or some other metric) between 
base and sensitivity biomass 
for each year in the model time series. It will take a directory input pointing 
to a particular sensitivity run (e.g. `data_outputs/sensitivity_1`) and returns
some summary statistics (such as mean percent error, mean absolute percent error, 
mins and maxs, etc.). As a side effect, this function will write the error metrics 
to a file (e.g. `data_outputs/sensitivity_1/sensitivity_comparison.csv`).

### Module 6. `plot_sensitivity_comparison.py` 

#### plots percent error between base and sensitivity biomass by year

This module will create a function which creates a time series plot of the error 
metrics written by `plot_sensitivity_comparison()`. It will take a directory input
(e.g. `data_outputs/sensitivity_1`), read `sensitivity_comparison.csv` in that 
directory, and either return a plot or save it to a file. It should probably test 
if the directory input and file exist and, if not, call `sensitivity_comparison()`.


## List of shiny app components to implement

### User interface

1. [ ] Natural mortality slider bar input
2. [ ] Natural mortality number box input
3. [ ] Other buttons/sliders/etc. for other parameters as place holders
4. [ ] Boxes to report error metric summary statistics
5. [ ] Box to show sensitivity/base biomass time series plot
6. [ ] Box to show error time series plot
7. [ ] Tabs to create multiple sensitivity runs

### Server

1. [ ] Take natural mortality input and call `modify_mortality()` then call `run_basa()`
2. [ ] Call `read_biomass_sensitivity()` reactive to `run_basa()`
3. [ ] Call `plot_sensitivity()` reactive to `read_biomass_sensitivity()` 
4. [ ] Call `sensitivity_comparison()` reactive to `read_biomass_sensitivity()` 
5. [ ] Call `plot_sensitivity_comparison()` reactive to `sensitivity_comparison()`
6. [ ] On creation of new tabs, new directories for sensitivity runs are created
       when `run_basa()` is called.
