## instructions to install BASA and its dependencies and execute the model

1. clone the repo

bash:

```
git clone git@github.com:cl-roberts/pws-herring-basa.git
cd pws-herring-basa
```

2. install `pwsHerringBasa`. This is an R package that I created which provides functions for compiling and executing the model, as well as analyzing its outputs.


bash:

```
R CMD INSTALL pwsHerringBasa_0.1.0.tar.gz
```

3. install other R package requirements.


4. Execute the model. This step should take a few minutes.

bash:

```
Rscript run_basa.r
```


5. Run analysis code. This should create some plots in `figures` and `.csv` files in `data_outputs`.

bash:

```
cd plotting
Rscript plot_management_outputs.r
Rscript plot_survey_fits.r
Rscript plot_age_compositions.r
Rscript plot_retrospective.r
```