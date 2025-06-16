This repository contains code from the paper "Multi-source analyses of average treatment effects with failure time outcomes", co-authored by Wen et al. 
The files contain codes to generate data that can be used in simulation studies, and TMLE methods to estimate causal estimands described in the main paper. The files in this repository include:

1. datagen.R: contains the function to generate data sets
2. tmle_withA6.R (tmle.R): codes to produce the parameter estimates from the TMLE methods described in the paper
(Note: tmle_withA6.R contains code for TMLE under Assumption 6, which assumes data source–outcome independence. tmle.R contains code for TMLE without this assumption)
3. mydata.csv: a sample dataset intended for use with the outlined analytical methods
