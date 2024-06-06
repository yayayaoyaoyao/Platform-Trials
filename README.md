# Platform-Trials

This repository provides some code for simulating response-adaptive randomisation with a Bayesian paradigm in platform trials to simultaneously test multiple interventions against a reference intervention within each subgroup in a factorial design. 

The following descriptions are from REMAP-CAP Team (2019a), REMAP-CAP Team (2019b), Mahar
et al. (2023), Tong et al. (2022) and REMAP-CAP Team (2020).

The simulation model can be simplified as
$$ log \left( \frac{p}{1-p}\right) = \alpha_{shock} + \beta_A + \beta_S + \beta_M + \delta_{shock,S2} + \delta_{shock,S3}  $$

with $p$ stands for the mortality rate, $\alpha_{shock}$ represents two strata (shock and non-shock), $A$ refers to the
antibiotic domain with 5 interventions, $S$ is the corticosteroid domain with 3 interventions, $M$ defines the
macrolide duration domain with 2 interventions, and $\delta$ is the interaction term between $S2$, $S3$ and shock
status. In this model, $A_{1−5}$ are nested interventions, which enable borrowing of estimated information
across interventions within the domain. So hierarchical priors are applied with $\beta_{A_{1−4}} \sim N(\mu, \tau^2)$, where
$\mu \sim N(0, 10^2)$ and $\tau^2 \sim IG(0.125, 0.00281)$. The effect of $\beta_{A_{5}}$ is set to 0 for identifiable purpose. 
Meanwhile, $S_{1−3}$ and $M_{1−2}$ are non-nested interventions, independent priors are used with $\beta_{S_{2,3}} \sim N(0, 10^2)$and
$\beta_{M_2} \sim N(0, 10^2)$, in which $S_1$ and $M_1$ are set to 0. The interaction between S and shock is modeled as
$\delta_{shock,S_2}, \delta_{shock,S_3} \sim N(0, 0.15^2)$.


JAGS (Just Another Gibbs Sampler) (Plummer et al., 2003; Plummer (2019), R implementation of JAGS) is a software library and a programming language from
which can be built Markov chain Monte Carlo (MCMC) samplers for Bayesian statistical models, but are slow with the complexity of model 
increasing and the number of participants increasing. We contribute a maximum likelihood approximation to posterior mean
to simulate a multi-factorial platform design with response-adaptive randomisation, with higher speed compared to JAGS. Such
findings could improve computation speed with acceptable accuracy, and it could provide valuable
insights for future research in platform trials. 

RAR_with_a_Bayesian_Paradigm_JAGS.R is the R code file using JAGS.

RAR_with_a_Bayesian_Paradigm_MLE.R is the R code file using maximum likelihood approximation.

functions.R is the R code file with required functions to run RAR_with_a_Bayesian_Paradigm_JAGS.R and RAR_with_a_Bayesian_Paradigm_MLE.R.

# References

REMAP-CAP Team (2019a). Appendix to Core Protocol: STATISTICAL ANALYSIS APPENDIX REMAP-CAP: Randomized, Embedded, Multifactorial Adaptive Platform trial for Community-Acquired Pneumonia. Tech. rep. URL: https://static1.squarespace.com/static/5cde3c7d9a
69340001d79ffe/t/5e82b68a143343364b38c10c/1585624721568/REMAP-CAP\+Statistical\+Analysis\+Appendix\+V3\+\-\+24\+August\+2019\_WM.pdf.

REMAP-CAP Team (2019b). Randomized, Embedded, Multifactorial Adaptive Platform trial for Community Acquired Pneumonia (REMAP-CAP): Core Protocol. Tech. rep. URL: https://static1.squarespace.com/static/5cde3c7d9a69340001d79ffe/t/5e82b62eda9c1c714065efc0/1585624649154/REMAP\-CAP\ +Core\+Protocol\+V3\+\-\+10\+July\+2019\_WM.pdf.

Mahar, R. K., McGlothlin, A., Dymock, M., Lee, T. C., Lewis, R. J., Lumley, T., Mora, J., Price, D. J.,
Saville, B. R., Snelling, T., et al. (2023). “A blueprint for a multi-disease, multi-domain Bayesian adaptive platform trial incorporating adult and paediatric subgroups: the Staphylococcus aureus Network
Adaptive Platform Trial”. Trials, 24(1), pp. 1–15.

Tong, S. Y., Mora, J., Bowen, A. C., Cheng, M. P., Daneman, N., Goodman, A. L., Heriot, G. S., Lee,
T. C., Lewis, R. J., Lye, D. C., et al. (2022). “The Staphylococcus aureus network adaptive platform
trial protocol: new tools for an old foe”. Clinical Infectious Diseases, 75(11), pp. 2027–2034.

REMAP-CAP Team (2020). Clinical Simulation Report: STATISTICAL ANALYSIS APPENDIX REMAPpCAP: Randomized, Embedded, Multifactorial Adaptive Platform trial for Community-Acquired Pneumonia. Tech. rep. URL: https : //static1.squarespace.com/static/5cde3c7d9a69340001d79ffe/t/ 5e8ed907493fb4203142d0fb/1586420012673/REMAP-CAP+SimulationReport-06+April+2020.pdf.

Plummer, M. et al. (2003). “JAGS: A program for analysis of Bayesian graphical models using Gibbs
sampling”. In: Proceedings of the 3rd International Workshop on Distributed Statistical Computing,
pp. 1–10.

Plummer, M. (2019). rjags: Bayesian Graphical Models using MCMC. R package version 4-10. 