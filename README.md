# Platform-Trials

This repository provides some code for simulating response-adaptive randomisation with a Bayesian paradigm in platform trials to simultaneously test
multiple interventions against a reference intervention within each subgroup in a factorial design. 

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

Plummer, M. et al. (2003). “JAGS: A program for analysis of Bayesian graphical models using Gibbs
sampling”. In: Proceedings of the 3rd International Workshop on Distributed Statistical Computing,
pp. 1–10.

Plummer, M. (2019). rjags: Bayesian Graphical Models using MCMC. R package version 4-10. 