Code for the leading-eigenvalue approximation to distribution of quadratic forms. Still needs to be packaged up.

fht.c: fast Hadamard transform
skat-big.R: Lanczos algorithm, using the 'svd' package
srht.R: algorithms based on sampling, both stochastic SVD and Lanczos algorithm with randomized trace estimator.

skat-sim2: example of simulations with the code.  Uses Gary Chen's MaCS simulator (https://github.com/gchen98/macs) to make the data

ToDo: 
  - move the computation of the test statistic and matrices into functions
  - automate the sparse-matrix construction
  - put it all in a package
  - make it work with seqMeta (?SKAT, ?RAREMETAL)
  
