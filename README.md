---
output: github_document
---


## A Case for Beta Regression in the Natural Sciences  
In preparation for submission to Ecosphere

* Authors:
    + [Emilie A. Geissinger](https://eageissinger.github.io/)
    + Celyn L.L. Khoo
    + [Isabella C. Richmond](https://github.com/icrichmond)
    + Sally J.M. Faulkner
    + [David C. Schneider](https://www.mun.ca/osc/dschneider/bio.php)
    

This repository contains the code and data accompanying the paper "A Case for Beta Regression in the Natural Sciences" (In Prep.). R scripts are organized in `scripts/`, raw data can be found in 'input/', figures used in the manuscript can be found in 'figures/', and model output tables can be found in 'output/'. Package dependencies include: `betareg`, `lmtest`, `DescTools`, `dplyr`.


    
## Abstract
Data in the natural sciences often take the form of percentages or proportions that are continuous and bounded by 0 and 1. Statistical analysis assuming a normal error structure can produce biased estimates when data are doubly bound. Beta regression uses an error structure appropriate for such data. We present examples from two fields of the natural sciences that produce continuous, bounded data: analytical chemistry and ecological elemental composition to compare the results from beta regression to normal error regression with arcsin transformation. We found that beta regression produced better model diagnostics, similar likelihood ratios and p-values, and substantially different parameter estimates. Our literature review showed that beta regression is used for a small portion of studies with data bounded at 0 and 1, with only a few examples in natural sciences. In our two case studies, normal error with arcsin transformation produced substantially different parameter estimates compared to beta regression. Beta regression for proportion data produced better model diagnostics and parameter estimates lying within the bounds of the data. The result is more accurate predictions and more reliable parameter estimates for use by other researchers. 
