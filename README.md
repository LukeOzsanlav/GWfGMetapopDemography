# Code for 'Differing drivers of decline within a metapopulation has implications for future conservation'

This repository holds the code for a publication that is submitted to *Ecology and Evolution* entitled; 'Differing drivers of reproductive output within metapopulations has conservation implications'

## _Authors_:

- Luke Ozsanlav-Harris <a itemprop="sameAs" content="https://orcid.org/0000-0003-3889-6722" href="https://orcid.org/0000-0003-3889-6722" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Geoff Hilton <a itemprop="sameAs" content="https://orcid.org/0000-0001-9062-3030" href="https://orcid.org/0000-0001-9062-3030" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
- Larry Griffin
- Lei Cao
- Mitch Weegman <a itemprop="sameAs" content="https://orcid.org/0000-0003-1633-0920" href="https://orcid.org/0000-0003-1633-0920" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
- Stuart Bearhop <a itemprop="sameAs" content="https://orcid.org/0000-0002-5864-0129" href="https://orcid.org/0000-0002-5864-0129" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

## Data Availability

Please note that the tracking data set used in the individual-level analysis has not been made available. However, the outputs of each script have been included so our analysis should be mostly repeatable. This means that `scripts 1-4` will not run by downloading the repository. The population-level analysis will run in full and the environmental data has been provided but the specific locations of nest sites removed.


## Manuscript Status

A pre-print of the Manuscript can be found [here](https://www.authorea.com/users/574634/articles/618299-differing-drivers-of-decline-within-a-metapopulation-has-implications-for-future-conservation?commit=5418ddc0fe5199421f4b27a6c2b70ca0a852796c)

MS submitted to *Journal of Animal Ecology* 19/12/2022

MS resubmitted to *Ecology and Evlution* 07/01/2023
Major revisions received 15/03/2023



## [Code description](Code)

- `1- Individual- Phenology calculator.R`: Calculates the departure and arrival dates into Greenland and Iceland for each tracked bird                                 
- `2- Individual- Classify incubations.R`: Calculates the incubation length for each tracked bird                     
- `3- Individual- Calculate nest locations.R`: From the incubation period this script calculates the putative nest location                 
- `4- Individual- Derive explanatory variables from Env-data.R`: Calculate the climatic conditions that each tracked bird experiences on the breeding grounds 
- `5- Individual- Create dataframe for models.R`: Create the data frame for individual-level models             
- `6- Individual- Breeding success & propensity models.R`: Run the binary GLMMs for breeding propensity and success    
- `7- Individual- Nest survival models.R`: Run cox ph models for nest survival                     
- `8- Individual- Repeatability of Greenland Arrival date.R`: Use rptR package to calculate repeatability of Greenland arrival dates
- `9- Population- Create pop-level sampling sites.R`: Choose the sites where sub-population climactic conditions are calculated
- `10- Population- Climatic effects on pop-productivity.R`: Model how the temperature and precipitation correlate with population-level productivity    
- `11- Create population trend graphs.R`: Creates figure 1 in the MS
- `Helper functions`: Additional functions to call for modelling


## Abstract

1.	To understand metapopulation declines, researchers generally ascribe demographic drivers in a single or few sub-populations and presume they are widely representative. With this information, practitioners implement blanket conservation measures across metapopulations to reverse declines by targeting the proximate mechanism(s) for overall population change. However, such approaches may not be appropriate in circumstances where sub-populations are spatiotemporally segregated and exposed to different environmental variation. 
2.	The Greenland White-fronted Goose Anser albifrons flavirostris is an Arctic-nesting migrant that largely comprises two sub-populations, which are delineated by northerly and southerly breeding areas in west Greenland. The metapopulation has declined since 1999 but this trend is only mirrored in one sub-population and the causes of this disparity are unclear. Here we compare the drivers and trends of productivity in both sub-populations using population- and individual-level analysis. 
3.	We examined how broad-scale temperature and precipitation influenced population-level reproductive success and whether there was a change in the relationship when metapopulation decline commenced. In addition we used biologging devices to reconstruct incubation events and modelled how phenology and environmental conditions influenced individual-level nest survival. 
4.	Correlations between reproductive success and temperature/precipitation on the breeding grounds have weakened for both sub-populations. This has resulted in declining reproductive success for the more northerly, but not the southerly breeding sub-population. At the individual-level this appears to be driven by differences in nest failure. Nest survival of the northerly breeding sub-population also increased with earlier breeding ground arrival and less precipitation during incubation, while no factors examined were important for the southerly breeding sub-population. This suggests reproductive success in the two sub-populations is now limited by different factors and that survival and/or emigration may be a more important driver of population trend for the southerly breeding sub-population.
5.	Demographic rates and their environmental drivers differ between the sub-populations examined here and consequently we encourage further decomposition of demography within metapopulations. This is important for conservation practitioners to consider more broadly as bespoke conservation strategies, targeting different limiting factors, may be required for different sub-population. 
