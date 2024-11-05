
## User stories

### User type 1 - PhD Student

#### CL Roberts

CL is a PhD student who will be the primary user of the project. He will use the sensitivity tool to inform aspects of his dissertation. He needs clear outputs from the tool to inform how BASA may be improved. He has a high technical skill level.

### User type 2 - Faculty

#### Trevor Branch

Trevor Branch is CL's PhD supervisor. He will use the sensitivity tool to inform aspects of his dissertation.


### User type 3 - Ecosystem researchers in Prince William Sound

The third type of user of this tool includes miscellaneous researchers in Prince William Sound - zooplankton, marine mammal, sea bird scientists, etc. - who use the herring model for their own studies. This user base will be able to use the tool to inform how various model assumptions may affect conclusions from their own research. They will need the sensitivity tool to be extensible to parameters that are most relevant to their research (i.e., model assumptions about herring reproduction will be relevant to pink salmon studies because pink salmon predate larval herring). This group may or may not be technical, but they have expert-level knowledge in their respective domains.

### User type 4 - Board members of agency funding model development
### User type 5 - General public

## TODO: write specific user stories for these types of users.

### Technician users

CL in 1 year will be the primary person who services this tool.


## Use Cases

#### Explicit use case
- Exploring the analysis (both graphically and maybe analytically)
- Exporting the results (both data files and graphs)
- Being able to save multiple runs for comparison

#### Implied use case

- Iteratively using the tool on an annual basis.


## Components

#### Adjust mortality slider

Name: Adjust mortality slider
What it does: Takes in an adjusted value for the instantaneous natural mortality assumption and calls on a Python module
Inputs (with type information): Adjusted assumption value for instantaneous natural mortality as provided by user (float64)
Outputs (with type information): None, directly
Other components used: None
Side effects: TBD

#### Adjust mortality module

Name: Adjust mortality module
What it does: Runs a Python module to produce a sensitivity-adjusted median biomass dataset
Inputs (with type information): Adjusted assumption value for instantaneous natural mortality as provided by user (float64)
Outputs (with type information): A sensitivity-adjusted median biomass dataset (dataframe)
Other components used: None
Side effects: May save over an existing sensitivity-adjusted dataset if one has been previously created

#### Residual time plot

Name: Create a residual time plot
What it does: Runs a Python module to create a matplotlib plot of the residuals between original and sensitivity-adjusted median biomass values.
Inputs (with type information): Original median biomass values (dataframe), adjusted assumption value for instantaneous natural mortality as provided by user (float64)
Outputs (with type information): A matplotlib plot 
Other components used: Adjust mortality slide, adjust mortality module
Side effects: May save over an existing sensitivity-adjusted dataset if one has been previously created

#### Calculate statistics

Name: Produce measures of difference
What it does: Runs a Python module to produce relevant measures of difference, such as root mean squared error, for the residuals between original and sensitivity-adjusted median biomass values
Inputs (with type information): Original median biomass values (dataframe), adjusted assumption value for instantaneous natural mortality as provided by user (float64)
Outputs (with type information): A list of measures of difference (dataframe)
Other components used: Adjust mortality slider, adjust mortality module
Side effects: May save over an existing sensitivity-adjusted dataset if one has been previously created

#### Exporting results

Name: Exporting/downloading results
What it does: Allows the user to save previous runs, including the calculated statistics and associated graphs/figures, to their computer
Inputs (with type information): internal sensitivity run data
Outputs (with type information): CSV files with the results that downloads to the user’s computer
Components used
Side effects: TBD

#### Cache sensitivity runs

Name: Saving results in a separate tab
What it does: Allows for comparing different sensitivity results WITHIN the webapp itself (split screens?)
Inputs (with type information): multiple sensitivity runs
Outputs (with type information): a tab for each sensitivity run that the user can click between 
Side effects:
Everytime you run the model, the package will generate a lot of plots (need to flesh out where these are saved and NOT saved (i.e. ignoring them on the git repo))

