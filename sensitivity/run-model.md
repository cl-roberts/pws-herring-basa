## instructions to install BASA and its dependencies and execute the model

1. clone the repo and switch to the `sensitivity` branch. Note that we will do all of our work in this branch

bash:

```
git clone git@github.com:cl-roberts/pws-herring-basa.git
cd pws-herring-basa
git switch sensitivity
```

2. install `pwsHerringBasa`. This is an R package that I created which provides functions for compiling and executing the model, as well as analyzing its outputs.


bash:

```
R CMD INSTALL pwsHerringBasa_0.1.0.tar.gz
```

Note that this step may prompt you to acquire some packages or linux libraries. For example, I needed to run `sudo apt install libcurl4-ssl-dev` to get this to install correctly. 

3. install other R package requirements. You should have already installed a vector of packages that I sent you in an `.rds` file, here are some others that were not included there.

R:

```
install.packages(c("TMB", "tmbstan"))
```

4. Execute the model. This step should take a few minutes.

bash:

```
Rscript run_basa_tmb.r
```


5. Run analysis code. This should create some plots in `figures/tmb` and `.csv` files in `data_outputs`.

bash:

```
cd plotting
Rscript plot_management_outputs.r
```
