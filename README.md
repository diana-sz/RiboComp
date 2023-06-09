# Ribosome composition

Scripts to study the effect of ribosome composition on growth rate using resource balance analysis (RBA) and elementary growth vectors


## code/

### Simulations
* general.py -- models, functions and classes used by other scripts 
* run_all_conditions.py -- RBA simulations that find max. growth rate at each ribosome composition
* fluxes_vs_growth_rate.py -- vary growth rate at fixed xrP = 36%; save RNAP fluxes
                              with or without accumulation of R or rRNA
* kostinski_reproduction.py -- RBA with fixed allocation of R and RNAP
                               parameters from Kostinski & Reuveni 2020
* fit_kdegmax.py -- fit R degradation rate so that the optimal composition is xrP=36%

### Plots 
* plot_xrp_vs_mu.R -- plot growth rate vs. ribosome composition
* plot_allocation.R -- plot ribosome allocation vs. ribosome composition
* plot_RNAP_fluxes.R -- plot RNAP fluxes vs. growth rate at fixed xrP=36%
* plot_RP_ratio.R -- plot RNA:protein ratio at different growth rates
* plot_vdeg_ratio.R -- plot ratio of RNA degradation flux to RNAP flux


## data/
* Outputs of the simulations
* * RBA_deg_[].csv -- output of run_all_conditions.py
* * fluxes_x0.36.csv -- output of fluxes_vs_growth_rate.py
* * Kostinski_reproduction_growth_rates.csv -- output of kostinski_reproduction.py
* parameters.csv - parameters for simulations in run_all_conditions.py
* gausing_RNA_deg.csv - fraction of degraded RNA at different growth rates from Gausing 1977 (extracted with WebPlotDigitizer)
* fluxes_bremer.csv -- fluxes from Bremer 1996, converted to mmol/gh with fluxes_vs_growth_rate.py


## plots/
Plots generated by R scripts in the section 'Plots'
