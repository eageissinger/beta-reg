## A Case for Beta Regression in the Natural Sciences  
In preparation for publication in Ecosphere

* Authors:
    + [Emilie A. Geissinger](https://eageissinger.github.io/)
    + Celyn L.L. Khoo
    + [Isabella C. Richmond](https://github.com/icrichmond)
    + Sally J.M. Faulkner
    + [David C. Schneider](https://www.mun.ca/osc/dschneider/bio.php)
    

This repository contains the code and data accompanying the paper "A Case for Beta Regression in the Natural Sciences" (In Prep.). R scripts are organized in `scripts/`, raw data can be found in `input/`, figures used in the manuscript can be found in `figures/`, and model output tables can be found in `output/`.  

Package dependencies include: `betareg`, `lmtest`, `DescTools`, `dplyr`.


    
## Abstract
Data in the natural sciences are often in the form of percentages or proportions that are continuous and bounded by 0 and 1. Statistical analysis assuming a normal error structure can produce biased and incorrect estimates when data are doubly bound. Beta regression uses an error structure appropriate for such data. We conducted a literature review of percent and proportion data from 2004 to 2020 to determine the types of analyses used for (0,1) bounded data. Our literature review showed that before 2012, angular transformations accounted for 93% of analyses of proportion or percent data. After 2012, angular transformation accounted for 52% of analyses and beta regression accounted for 14% of analyses. We compared a linear model with angular transformation with beta regression using data from two fields of the natural sciences that produce continuous, bounded data: biogeochemistry and ecological elemental composition. We found little difference in model diagnostics, likelihood ratios, and p-values between the two models. However, we found different coefficient estimates from the back-calculated beta regression and angular transformation models. Beta regression provides reliable parameter estimates in natural science studies where effect sizes are considered as important as hypothesis testing.
