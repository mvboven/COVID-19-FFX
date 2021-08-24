Data and scripts of a SARS-Cov-2 household transmission study as described in DFM Reukers, M van Boven, A Meijer, N Rots, C Reusken, I Roof, AB van Gageldonk-Lafeber, W van der Hoek, S van den Hof, High infection secondary attack rates of SARS-CoV-2 in Dutch households revealed by dense sampling, Clinical Infectious Diseases (2021), https://doi.org/10.1093/cid/ciab237.

The data directory contains the data of the household study. Date of birth of the participants and date of first symptoms have been removed from the data for privacy reasons. The data directoy also contains a codebook for the data (in Dutch).

The scripts directory contains an R script with which the analyses in the paper can be reproduced. The first part contains code for the logistic regressions using Generalizsed Estimating Equations. The second part of the R script contains code for analysing the data using multitype final size equations of SEIR type transmission models. The core model is implemented in Stan, and also available in the scripts directory.

Enquiries about the data should directed to Daphne Reukers (daphne.reukers@rivm.nl) and questions about the statistical analyses to Michiel van Boven (michiel.van.boven@rivm.nl).
