
## User stories

### User type 1 - PhD Student

#### CL Roberts

CL is a PhD student who will be the primary user of the project. He will use the sensitivity tool to inform aspects of his dissertation. He needs clear outputs from the tool to inform how BASA may be improved. He has a high technical skill level.

### User type 2 - Faculty

#### Trevor Branch

Trevor Branch is CL's PhD supervisor. He will use the sensitivity tool to inform aspects of his dissertation.


### User type 3 - Ecosystem researchers in Prince William Sound

The third type of user of this tool includes miscellaneous researchers in Prince William Sound - zooplankton, marine mammal, sea bird scientists, etc. - who use the herring model for their own studies. This user base will be able to use the tool to inform how various model assumptions may affect conclusions from their own research. They will need the sensitivity tool to be extensible to parameters that are most relevant to their research (i.e., model assumptions about herring reproduction will be relevant to pink salmon studies because pink salmon predate larval herring). This group may or may not be technical, but they have expert-level knowledge in their respective domains.

Mary studies humpback whales in Prince William Sound, where herring are a major food source for humpback whales. She wants to learn more about herring population for her research on humpback whales, and to assess how sensitive the model is to error to see how the food source for humpback whales may change. Mary wants to be able to see the biomass of herring without having to sort through too much data on her own. Mary has high technical skills but wants to be able to use the model with minimal effort.

James studies plankton in Prince William Sound, which is a major food source for herring. He wants to know more about the BASA models sensitivity for results regarding herring biomass to understand how these changes might affect plankton dynamics. James has intermediate to high technical skills but prefers a more intuitive interface to enable him to explore the data easily. 

Kiersten studies pink salmon in Prince William Sound, an ecological competitor of herring. She is particularly interested in how assumptions about herring influence biomass predictions to inform her understanding of competition dynamics. Kiersten has advanced technical skills and wants a flexible tool so she can run tailored sensitivity analyses. She also wants to be able to compare her results.

### User type 4 - Board members of agency funding model development

Sylvia serves as a board member for the agency funding the model's development as a member from the general public. As a board member, she wants to know how accurate the model is so she can best represent the best interests of the community. She has limited technical knowledge and values clear, accessible explanations and figures for the data.

Nancy serves as a board member for the agency funding the model's development. She wants to understand how accurate the model is to ensure the interests of commercial fishing are best represented. She wants a tool that will allow her to see how changes in variables affect the predictions with little effort on her part. Nancy has low to intermediate technical skills.

Bob serves as a board member for the agency funding the model's development. He wants to know how accurate the model is so he can best represent sustenance fishers. Bob wants a tool that simplifies the complex analyses down to laymans terms and clearly shows the sensitivity. Bob has low technical skills. 

### User type 5 - General public

Thomas is a fisherman in the Prince William Sound area. He is curious about how fishing populations would change 
if certain conditions were met (e.g. birth rates were higher or lower, if he made the decision to throw back older fish). 
He is very familiar with fishing terms, but not familiar at all with statistics terms (e.g. standard deviation, residuals). 
He has average tech literacy and he prefers a straight forward term as opposed to something with all the bells and whistles. 

Kerry is a 10th grade student living near the Prince William Sound area. He stumbled upon this as a random web app 
and is just curious to play around with it. He’s not familiar with fishing terms or statistics terms, 
but has a fairly high tech literacy. He won’t spend a lot of time using the web-app but is curious to 
play with it for a few minutes.

Lily is a college student studying ecology. She is working on a class project on overfishing. 
She would like to use the web-app to generate some hypothetical scenarios to include in her class project. 
She has high tech literacy, is fairly familiar with statistics terms, and a little bit familiar with fishing terms. 
It’s important that she can compare different sensitivity runs. 


## TODO: write specific user stories for these types of users.

### Technician users

CL in 1 year will be the primary person who services this tool.


## Use Cases

#### Explicit use case
- To generate the adjusted model output from the base data and a user-specified mortality value
- To generate a data series of residuals between base data and adjusted data
- To explore the difference between base model output and adjusted model output graphically, using a plot of time on the x axis and residuals between base data and adjusted data on the y axis (Plot A)
- To explore the difference between base model output and adjusted model output graphically, using a plot of time on the x axis and percent difference between base data and adjusted data on the y axis (Plot B)
- To visualize Plot A and Plot B side-by-side
- To generate analysis values from the base data and the adjusted model output, specifically root mean squared error (RMSE) and mean absolute error (MAE) to evaluate residuals
- To generate analysis values from the base data and the adjusted model output, specifically R^2 to evaluate how much of the adjusted model output variation can be explained by the original data, and therefore understand how much the adjusted model dataset is affected by the changing mortality value
- To export the adjusted model output and residual data
- To export image files for Plot A and Plot B
- To save a previous model run inside the app



#### Implied use case

- Run the model every year using the annually updated dataset


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

