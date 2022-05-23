
# Practicum project on microsimulation modeling for alcohol policy in Ontario

Author: Sze Hang (Hana) Fu (hana.fu@mail.utoronto.ca)
Supervisor: Lennon Li (Lennon.Li@oahpp.ca) 
Co-supervisors: Brendan Smith (brendan.smith@oahpp.ca) and Erin Hobin (erin.hobin@oahpp.ca)

*Last Updated: May 23, 2022*

## About

This GitHub repository contains R scripts and other materials for the practicum placement at Public Health Ontario (PHO), during October 2021 to April 2022. The focus on the practicum placement project is on microsimulation modeling for alcohol policy in Ontario.

## Alcohol and microsimulation

Alcohol policy, from how alcohol is sold, taxed and priced, is continually changing in Ontario. The interest is to build an individual-level discrete-event microsimulation model, representative of the population of Ontario, to study the impact of changes to alcohol policies on population and social inequities in alcohol-attributable harms. 

Since we are at the initial stage of the project, the aim of the practicum placement is to set up basic structure for the alcohol policy microsimulation model to understand alcohol harms for the population in Ontario from 2013-2043. This initial model is used to predict alcohol-related mortality without change in policy effect, while consider the differences in age and sex for the Ontario population.

More details on this basic microsimulation structures can be found in the PowerPoint presentations in the [presentations folder](presentations). Also see the [Reference section](#reference) for reference papers on microsimulation modeling. 

## Folder structures

| Folder | Description | 
|------|-------|
| [Analysis](Analysis) | Contains R scripts for loading Canadian Community Health Survey (CCHS), organizing data for performing microsimulations, and performing the analysis. | 
| [Microsim_functions](Microsim_functions) | This folder contains functions in R scripts for performing microsimulations in R. The codes are based on the tutorial paper by Krijkamp EM et al. 2018 [https://doi.org/10.1177/0272989X18754513](https://doi.org/10.1177/0272989X18754513) and their R codes can be found in [https://github.com/DARTH-git/Microsimulation-tutorial](https://github.com/DARTH-git/Microsimulation-tutorial).  | 
| [Summaries](Summaries) | With R scripts for creating summary reports. | 
| [Presentations](Presentations) | Contains PowerPoint presentations for presentation in class and to PHO team. | 
| [openMpp](openMpp) | Contains R scripts for setting up analysis using openM++. openM++ is an open source microsimulation platform, based on the Modgen software ([https://www.statcan.gc.ca/en/microsimulation/modgen/modgen](https://www.statcan.gc.ca/en/microsimulation/modgen/modgen)). See [https://github.com/openmpp/R](https://github.com/openmpp/R) for more detail.| 

## References

Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. Microsimulation modeling for health decision sciences using R: A tutorial. Med Decis Making. 2018;38(3):400–22.


```bibtex
@article{
    jha2014,
    title={Reliable direct measurement of causes of death in low- and middle-income countries},
    author={Jha, Prabhat},
    journal={BMC Medicine},
    volume={12},
    number={1},
    pages={1--10},
    year={2014}
    DOI={10.1186/1741-7015-12-19},
    publisher={Springer}
}
```
