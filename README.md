# Empirical-Results-of-MaxMI-for-RRB-Identification

This project's main purpose is to provide R code that can be used to reproduce all empirical results in 'Mutual Information as a Tool for Optimal Classification: Application to Identifying Rapid-Responding Behaviour' by Santeri Holopainen, Jari Metsämuuronen, Mikko-Jussi Laakso, and Janne Kujala.

This project is maintained by Santeri Holopainen (email: sjholo@utu.fi), Turku Research Institute for Learning Analytics, University of Turku, Finland.

The R code are divided into three sets of scripts:
1. R code for reproducing the results of applying and analysing the MaxMI methods in the PISA 2022 data set:
    1. user-defined functions (1.0_pisa2022_functions.R)
    2. data management (1.1_pisa2022_data_management.R)
    3. results (1.2_pisa2022_results.R)
2. R code for reproducing the results of examining the behaviour of MaxMI-22 under certain conditions according to the ILC-IRT framework at the population level:
    1. user-defined functions (2.0_ILC-IRT_population_level_functions.R)
    2. results (2.1_ILC-IRT_population_level_results.R)
3. R code for reproducing the results of investigating the usability of MaxMI-22 under certain conditions according to the ILC-IRT framework at the realised level:
    1. user-defined functions (3.0_ILC-IRT_realised_level_functions.R)
    2. results (3.1_ILC-IRT_realised_level_results.R)

Some notes:
- The R code should be ready to run as is, provided that you have downloaded the PISA 2022 data as well as installed the necessary R packages.
- The scripts for the user-defined functions are run with the 'source' function from the other scripts.
- The PISA 2022 data is available at the OECD website: https://www.oecd.org/en/data/datasets/pisa-2022-database.html
- When running the code for reproducing the PISA 2022 results, one should start with the data management script.
