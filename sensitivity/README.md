# A GUI for Sensitivity Analysis for the Prince William Sound Herring Bayesian Age-Structured Stock Assessment Model (BASA)

This is a GUI for **easy** and **visually appealing** sensitivity analysis in the 
BASA model, which estimates population parameters of a herring stock in Prince 
William Sound in Alaska. The sensitivity analysis explores how the BASA model 
time series estimates for biomass change if the assumed instantaneous rate of
natural mortality is modified. The shiny app calls the `pwsHerringBasa` R package
and `sensitivity` python package that are provided in this repository

## Setup

0. Ensure that Python has been downloaded through Anaconda or Miniconda and added 
   to your system path as relevant. 

1. Clone the `sensitivity` branch of the `pws-herring-basa` repository and work
   from the repository root directory

```
git clone -b sensitivity --single-branch https://github.com/cl-roberts/pws-herring-basa.git

cd pws-herring-basa
```

2. Obtain python dependencies

```
conda install pandas numpy matplotlib pytest
conda install -c conda-forge shiny
```

3. `pip` install the `sensitivity` package. This contains modules that
   interface with the shiny app for running the sensitivity analysis

```
cd sensitivity
pip install --no-deps -e .
cd ..
```

4. Obtain R dependencies

```
Rscript -e 'install.packages(c("ggplot2", "dplyr", "snow", "data.table", "adnuts", "snowfall", "rstan", "r4ss"), repos="https://cloud.r-project.org")'
```

5. Install the `pwsHerringBasa` R package from the `tar.gz` file provided in the 
   repository root directory. This package provides functions to execute BASA and
   clean its outputs.

```
R CMD INSTALL pwsHerringBasa_0.1.0.tar.gz
```

6. Run the shiny app

```
shiny run --reload --launch-browser sensitivity/sensitivity-app/app.py 
```

If you run into any issues, the first step recommended is to pip uninstall any 
of the relevant packages and reinstall them through conda.

**Note:** BASA is written in ADMB, a specialized C++ wrapper for fitting highly 
parameterized fisheries stock assessment models. The ADMB model is provided as 
a `.tpl` file, and ADMB compiles it to an executable. I have removed the compilation
step from the sensitivity analysis and, instead, provided the model executable 
file `sensitivity/model/PWS_ASA.exe` to remove ADMB as a dependency to running 
the sensitivity analysis shiny app. 


## Contributions

CL Roberts:
- Team lead for overall project
- Designed package and module structures
- Designed and coded Shiny app

Veronica Lee:
- Wrote explicit and implied use cases
- Wrote Module 5: sensitivity_comparison.py
- Wrote test_sensitivity_comparison.py
- Wrote Module 6: plot_sensitivity_comparison.py
- Wrote test_plot_sensitivity_comparison.py
- Did pylint code review for Modules 5 and 6
- Created powerpoint slides for final presentation

Mindy Dai:
- Wrote README
- Wrote modify_mortality.py and associated tests
- Assisted with pylint code review

Trinity Hinshaw:
- Wrote user stories
- Wrote Module 4: plot_sensitivity.py
- Wrote test_plot_sensitivity.py
- Did pylint code review for Module 4 and associated tests file
