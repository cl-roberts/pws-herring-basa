# A GUI for Sensitivity Analysis for the Prince William Sound Herring Bayesian Age-Structured Stock Assessment Model (BASA)

This is a GUI for **easy** and **visually appealing** sensitivity analysis in the 
BASA model, which estimates population parameters of a herring stock in Prince 
William Sound in Alaska. The sensitivity analysis explores how the BASA model 
time series estimates for biomass change if the assumed instantaneous rate of
natural mortality is modified. The shiny app calls the `pwsHerringBasa` R package
and `sensitivity` python package that are provided in this repository

## Future Directions

#### Expanding on the test modules 
- Implement tests for the Shiny app using the `shinytest` and `runTests()` functionality
- In the short term, we would focus on mainly Shiny app unit tests to make sure our basic features were working
- In the long term, we would probably work on [snapshot-based tests](https://shiny.posit.co/r/articles/improve/testing-overview/) to help diagnosis our issues on why the app is running on certain machines but not others
- BASA is a probabilistic model so our outputs have an inherent degree of randomness, as a result we found it difficult to implement one-shot tests for our `sensitivity_comparison` module. If we had more time, we would have wanted to figured out what would be appropriate tolerance ranges.
#### Adding additional sensitivity parameters
- Our original intention was to be able to change multiple parameters simultaneously so in addition to "natural mortality rate", we also wanted to change other parameters such as "birth rate", "disease seroprevalence", "age distribution", among others. 
- Due to time constraints, we were only able to implement "natural mortality rate", but we tried to write the code in a modular fashion that would allow us to add more parameters following the same pattern that we implemented "natural mortality rate".
- In future iterations, we would like the Shiny app to have more sliders so it can perform more nuanced sensitivity analyses.
#### Implement multiple tabs functionality
- It is very common to want to compare te results of different sensitivity analyses- in fact, one could argue that it's one of the primary purpose of sensitivity analyses. In our current iteration, users can see the result of one sensitivity analysis, but would need to open a new session if they wanted to run a new analysis while also viewing the initial one. 
- In future iterations, we would like to add a multiple tab function so that users can view and compare different sensitivty runs side by side visually. 
- Shiny allows for the generation of multiple tabs and panels fairly easily, though it may be more difficult to implement screen splitting so that the user can view plots side by side. 
#### Implement a download feature
- In future iterations, we would like to add a button so that users can easily download their current sensitivity run outputs and associated plots onto their local computer.
#### Speed Up Loading Time
- Currently, it takes roughly 1-2 minutes for our app to generate output for each input change. This, unfortunately, is primarily due to limitations from the BASA model that our app is based off of. Since it's a simulation based model, it needs to run and analyze a large number of simulations, which is computationally intensive. 
- However, we understand that this long loading time is a significant barrier for anyone outside of our primary user (who needs to use it for work and research) to actually want to use it. 
- While we anticipate that this is the hardest "Next Step" to tackle, our team discussed several possible solutions such as decreasing the number of simulations. This would decrease the precision of results, but would speed up the loading time. 
- In the very long term, rewriting the base package in a faster language such as Julia might be the best solution, but that is beyond the scope of this particular subpackage.


## Setup

0. Ensure that Python (>=3.10) been downloaded through Anaconda or Miniconda and added 
   to your system path as relevant. Also ensure that [R](https://cran.rstudio.com/) 
   (>=4.4.0) is installed and `Rscript.exe` is added to your system path.

```
PATH="$PATH:~/path/to/R/bin/"
```

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
- Updated README with Next Steps
- Wrote modify_mortality.py and associated tests
- Added assorted additional edge tests
- Assisted with pylint code review

Trinity Hinshaw:
- Wrote user stories
- Wrote Module 4: plot_sensitivity.py
- Wrote test_plot_sensitivity.py
- Did pylint code review for Module 4 and associated tests file
