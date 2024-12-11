# A GUI for Sensitivity Analysis for the Prince William Sound Herring Bayesian Age-Structured Stock Assessment Model (BASA)

This is a GUI for **easy** and **visually appealing** sensitivity analaysis in the BASA model, which measures level of herring stock in the Prince William Sound in Alaska. It is built on top of the BASA package and generates a Shiny app.

## Setup

Ensure that Python has been downloaded through Anaconda or Miniconda and added to your path as relevant. 

Clone the sensitivity branch using the following text in the command line: git clone -b sensitivity --single-branch https://github.com/cl-roberts/pws-herring-basa.git

Change directory into pws-herring-basa

If needed, run the following commands in the command line: conda install pandas, conda install numpy, conda install matplotlib, conda install pytest

Install Shiny for Python in the command line: conda install -c conda-forge shiny

Change directory into pws-herring-basa/sensitivity

Install the package 'sensitivity' in the command line. Include the period at the end of this text: pip install --no-deps -e .

Change directory back to pws-herring-basa

Load the Shiny app from the command line: shiny run --reload --launch-browser sensitivity/sensitivity-app/app.py 

If you run into any issues, the first step recommended is to pip uninstall any of the relevant packages and reinstall them through conda.

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
