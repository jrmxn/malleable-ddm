Directory contents
========

- ssm_def_base.m
	- An example of how to inherit from ssm_def.m and extend basic functionality.
- ssm_def_conflict.m
    - A more detailed example of how to inherit from ssm_def.m and extend basic functionality to model the Simon effect.
- test_compare_likelihoods.m
	- Fits a model to the data ('testing.csv'), and compares likelihood estimation for different methods.
- test_compare_methods_pdf.m
    - Re-evaluates response time distribution from a model (testing_fit.mat)  with different methods (analytical solution, transition matrix approach, and brute force approach).
- test_ssm_run.m
    - An example of a model fit that uses ssm_def_base.m
- test_ssm_run_conflict.m
    - An example of a model fit that uses ssm_def_conflict.m
- testing.csv
    - Data in testing.csv is a cut down version (3 subjects) of data from OpenFMRI, accession:ds000101.
- testing_fit.mat
	- see test_compare_methods_pdf.m
