# cluster-se
approach to elegantly compute cluster robust standard errors in two-sample-2-stage-least-square using STATA and analytic expression of variance estimator. Adaption/extension of the code for hetereskedasticity robust standard errors as provided by

Robust inference for the Two-Sample 2SLS estimator,
Economics Letters,
Volume 146,
2016,
Pages 50-54,
ISSN 0165-1765,
https://doi.org/10.1016/j.econlet.2016.06.033


Synthetic data generated, Monte Carlo study showing equivalence of analytical estimator and bootstrap variance provided, description found in dissertation of mine.

BEWARE: There exists no proof, that the cluster robust derivation is valid, although it is intuitive, and produces plausible results using synthetic data, in fact, analytical and bootstrap standard error are the same with 10,000 simulation turns and a few specifications of paramaters tested.

2018 Stefan Etgeton, 2016 David Pacini/Frank Windmeijer
Creative Commons Licence:  CC BY-NC-SA (Attribution-NonCommercial-ShareAlike)
https://creativecommons.org/licenses/by-nc-sa/4.0/

