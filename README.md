# Ribosome composition

Scripts to study the effect of ribosome composition on growth rate using Elementary growth modes (EGMs) and elementary growth vectors (EGVs)


## code/

### Simulations
* 00_RBA.py -- standard RBA
* 01_RBA_reverse.py -- standard RBA with parameters that make RNA more expensive than protein
* 02_RNAPmax.py -- RBA with a limit on total RNA polymerase (RNAP)
* 03_RNAPmax_act.py -- RBA with a limit on total RNAP + limits on the maximum Ribosome (R) and RNAP activity
* 04_RNAPmax_arch.py -- RBA with a limit on total RNAP with archaeal parameters
* 05_fluxes_RNAPmax.py -- RBA with a limit on total RNAP; fixed xR = 36%; saves RNAP fluxes
* 06_fluxes_RNAPmax_noacc.py -- RBA with a limit on total RNAP, no accumulation of R or rRNA; fixed xR = 36%; saves RNAP fluxes
* 07_fluxes_RNAPmax_noacc_act.py -- RBA with a limit on total RNAP, no accumulation of R or rRNA; limits on the maximum R and RNAP activity; fixed xR = 36%; saves RNAP fluxes
* 08_PRL_reproduction.py -- RBA with fixed allocation of R and RNAP and maximum activities, parameters from Kostinski & Reuveni 2020

### Plots 
* plot_xR_vs_mu.R
* plot_RNAP_fluxes.R
* plot_PRL.py


## data/
Outputs of the simulations


## plots/
Plots generated by R scripts in the section 'Plots'
