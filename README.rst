************
Introduction
************

:Author: James R. McIntosh
:Contact: 
:Web site:
:Github:
:Mailing list:
:Copyright:
:License:
:Version: 0.0.1

Todo
====
# ddm_bru not currently used, and needs fixing + converting into a pdf
# maybe implement a chi-square minimisation, but meh.
# Could somehow do multiple subject MCMC to make it hierarchical... but might be a lot of work.
# At the moment ddm_def_instance returns handle to a random distribution as well as a single draw from the distribution. The single draw is needed but... it's a bit ugly that it's defined separately.


Purpose
=======
To make a flexible DDM fitting tool.
The idea here is to enable quick modification of the DDM rather than to optimise fitting methods for a pre-existing formulation.

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

Comparison to other packages
============================


Installation
============


How to cite
===========


Getting started
===============

