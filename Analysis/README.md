# Analysis

ontains R scripts for loading Canadian Community Health Survey (CCHS), organizing data for performing microsimulations, and performing the analysis

## Files

* [cchs_sim_data.Rmd](cchs_sim_data.Rmd): Codes for loading CCHS data, examining sample and bootstrap weights, create factor variables for socioeconomic and alcohol variables. Also include tabulation of counts by variables, implement survey weights using the [survey](survey) package, and dummy transition probabilities and relative ratios.
* [cchs_data_trans_prob_sim.Rmd](cchs_data_trans_prob_sim.Rmd): Include partial codes from [cchs_sim_data.Rmd](cchs_sim_data.Rmd), and addition codes on:
  * comparing unweighted counts and survey weighted counts
  * methods to combine the alcoholic and non-alcoholic paths
  * codes to run microsimulation models, sourcing codes from the [Microsim_functions](../Microsim_functions) folder
  * summarizing microsimulation results
