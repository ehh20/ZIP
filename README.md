
# ZIP

This repository contains some functions for fitting what might be called 'ZIP-INGARCH(1,q)' (zero-inflated Poisson Integer-valued GARCH) models (with exogeneous covariates) to a network count process.

Brief description of file contents:

-   `Functions_inputs_descriptions` contains a 'glossary' of the arguments used in the different functions.

-   `ZIP.R` has functions for estimating parameters for the univariate case and for forecasting and evaluating model fit.

-   `MZIP.R` does the above for a (specified) network. 3 different levels of degrees of freedom can be specified for fitting parameters, based on the network structure.

-   `MZIP_helper_functions.R` contains 'helper functions' that are used within the functions in `MZIP.R`. Functions in this file probably do not need to be called directly.

-   `PAR_HMM.R` is a basic implementation of the order (1,q) INGARCH model described in @fokianos2009. This is then extended to a 2-state hidden (homogeneous) markov model, where the other state emitting from a Poisson distribution with low rate.