# Generalized-Score-Functions-for-Causal-Discovery

Copyright (c) 2017-2018 Biwei Huang

Causal structure learning by greedy equivalence search with generalized score functions (which is applicable to mixed continuous and discrete data, data with Gaussian or non-Gaussian distributions, linear or nonlinear causal mechanisms, and variables with multi-dimensionalities.)

### IMPORTANT FUNCTIONS:

function [Record] = GES(X,score_type,maxP,parameters)

* INPUT:
  * X: Data with T*D dimensions
  * score_type: the score function you want to use
     *     score_type = 1: cross-validated likelihood
     *     score_type = 2: marginal likelihood
  * maxP: allowed maximum number of parents when searching the graph
  * parameters: when using CV likelihood, 
     *      parameters.kfold: k-fold cross validation
     *      parameters.lambda: regularization parameter

* OUTPUT:
  * Record.G: learned causal graph
  * Record.update1: each update (Insert operator) in the forward step
  * Record.update2: each update (Delete operator) in the backward step
  * Record.G_step1: learned graph at each step in the forward step
  * Record.G_step2: learned graph at each step in the backward step
  * Record.score: the score of the learned graph


### EXAMPLE:
see example1.m


### CITATION:
	
Biwei Huang ,Kun Zhang , Yizhu Lin, Bernhard Scholkopf, Clark Glymour. Generalized Score Functions for Causal Discovery. KDD, 2018.
