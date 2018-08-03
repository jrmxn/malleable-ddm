************
Introduction
************

:Author: James R. McIntosh
:Contact: j.mcintosh@columbia.edu
:Web site:
:Github:
:Mailing list:
:Copyright:
:License:
:Version: 0.0.5

Purpose
=======

To allow rapid prototyping of DDM extensions.
Code is a synthesis of tricks that I have been using to run drift-diffusion based models for a few years, as well as some taken from the HDDM toolbox.

Features
========

- Easily extensible framework:
	- Inherit from a single function (ddm_def.m)
	- Define all your additional model parameters
	- Specify the link between data (e.g. trial-trial eeg, visual stimuli) and model inputs
	- Write your likelihood function making use of model parameters and model inputs (examples provided)
- Because this is all specified directly, the model structure is very clear
- Model is specified as two binary numbers representing each model parameter. One number represents all the set parameters for this evaluation, the other number represents all the parameters over which to optimise.
e.g. if our full model is s,a,t,v,st then we could define the following submodel: model_id = [1,1,1,1,0], model_search = [0,1,1,1,0], which would initialise a search by setting st to 0 (unused), s to its default value (1) and randomly initialise and then optimise a, t and v.
	- This means it's easy to setup complex models and systematically search submodels, or freeze certain parts of the model while allowing some parameters to be free.
	- Simple to start optimisation of complex models based on simpler models
	- Easy to keep track of many candidate moels
- MCMC estimation available, and cost function minimisation uses a pattern search algorithm (easily modifiable)
- Written in Matlab (a feature for some!). Minimisation of cost function
- Fits are saved as matlab objects - makes it very easy to load previous fits, extend, re-fit, and keep track of the model fitting history.
- Three likelihood calculation methods are currently included as examples:
	- Core functions of HDDM (Navarro and Fuss analytical solutions) for speed - when initialising models or only using simple models
	- Brute force simulation with density estimation (using matlab's ksdensity) - easy to extend, but slow and prone to local minima due to stochastic nature (although pattern search generally does a good job)
	- Calculation decision-variable path with transition matrix iterative multiplication - easy to extend, quite fast with reasonable discretization (apart from for specific parameters, e.g. including noise in drift is a x10 slow down).

Notes on function
========

Without using an analytical solution we do not have the advantage of being able to generate likelihoods without generating the whole p(rt,accuracy|parameters,conditions).
For each trial we consequently create a structure of parameters and conditions e.g. parameters: p(ix_trial).s = 1;p(ix_trial).v = 1; condition:p.(ix_trial).conflict = 1;
Then this structure is stripped of duplicates, the likelihood density is estimated, and the relevant accuracies/rt are evaluated from it.

Weaknesses
========
It's not hierarchical.
Straying from the DDM means you lose the advantage of analytical solutions - i.e. things get slow(er)!

Features
========
For example, adding a conflict dependent bias is straight forwards:

1) inherit from ddm_def

2) modify the method ddm_cost_add_stim_dependencies so that the core ddm functions can make use of conflict
    (e.g. one line: p_mat.c = obj.data.stim_conflict;)
	
3) overload ddm_def_instance with the introduced bias parameter. Set its prior, default and bounds.

4) write the core ddm_pdf function (the current implemented version uses a transition matrix approach)

The structure is like this: the modelclass defines the space of all parameters you want to include in this specific version of the DDM.
Then there is an id_model and id_search.
id_model sets the specific subset of parameters that you want to be non-default (I hope not non-zero).
id_search sets the which of these parameters are then optimised (for example, you want noise in the model, but you don't want it to be a free parameter).
id_model and id_search are fundamentally binary strings that represent the inclusion of the model, and are converted to decimals and used to save/restore files.


Installation
============


How to cite
===========
TBD - please check back in a couple of months!

Getting started
===============


Todo
====
- Could somehow do multiple subject MCMC to make it hierarchical... but might be a lot of work.