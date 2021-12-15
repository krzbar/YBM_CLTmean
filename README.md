# YBM_CLTmean
These are the R scripts and numerical results accompanying "A Central Limit Theorem for branching Brownian motion".

The R setup for the manuscript was as follows: R version 3.6.1 (2019-09-12) Platform: x86_64-pc-linux-gnu (64-bit) Running under: openSUSE Leap 42.3

The exact output can depend on the random seed. However, in the script we have the option of rerunning the analyses as it was in the manuscript, i.e. the random seeds that were used to generate the results are saved, included and can be read in.

To replicate the results of the simulation study one needs to first in each of the folders (for i=30, 100, 1000, 10000, 20000) YBMlim_simulations/n_i unzip the zip archive RandomSeeds_YBMsimulstudy_n_i.zip. It is important that the directory PartialResults in the zip archive is placed in the directory n_i . After the random seeds are unpacked one runs source("torun_YBMlim_simul.R"). In order to replicate the results the logical variable b_use_random_seed has to be set to TRUE. If not, then the simulations will start from the provided by R .Random.seed. If furthermore b_save_random_seeds is TRUE then the random number generator seeds will be saved (hence a combination of FALSE and TRUE respectively will lead to overwriting whatever random number generator seeds are stored. The integer variable numcores is the number of CPU cores to be used in the parallalization and runnum is a (character) suffix that will be added to file and directory names if one wants to have multiple runs with distinct names that do not overwrite each other.
