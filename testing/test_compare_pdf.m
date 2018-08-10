clear;
addpath('..');
%%
% Compare likelihood evaluations for DDM without extensions
% analytical solution from NF/HDDM is most accurate
% transition matrix approach and brute force approach approximate
% analytical solution closely.
% With sufficient iterations, brute force approach and transition matrix
% approach are virtually identical - although they deviate slightly from
% analytical solution unless time resolution (dt) is sufficiently low.
%%
sr = load('testing_fit.mat');sr = sr.obj;%just load this as a template
p = sr.fit(end).p;
sr.s.dt = 5e-4;%for transition matrix and brute force approach
sr.s.nits = 25e3;%for brute force approach
sr.s.dx = 5e-3;%for transition matrix approach
%%
lt = [0:sr.s.dt:5];
p.t = 0.25;
p.s = 1;
p.c = 0;
p.v = 1.5;
p.st = 0.00;
p.sz = 0.00;
p.sv = 0.0;
p.z = 0.5;
p.a = 1.0;
%%
tic
[pdf_cw_ana,pdf_cr_ana,rt_ana] = sr.ddm_pdf_ana(p,lt);
toc
%%
tic
[pdf_cw_trm,pdf_cr_trm,rt_trm] = sr.ddm_pdf_trm(p,lt,sr.s.dx);
toc
%%
tic
[pdf_cw_bru,pdf_cr_bru,rt_bru,cdf_dow,cdf_ups] = sr.ddm_pdf_bru(p,lt,sr.s.nits);
toc
%%
close all;clf;hold on;
lw = 1;

ar = plot(rt_bru,pdf_cr_bru,'r','lineWidth',lw);
plot(rt_bru,pdf_cw_bru,'--r','lineWidth',lw);

ag = plot(rt_trm,pdf_cr_trm,'g','lineWidth',lw);
plot(rt_trm,pdf_cw_trm,'--g','lineWidth',lw);

ab = plot(rt_ana,pdf_cw_ana,'b','lineWidth',lw);
plot(rt_ana,pdf_cr_ana,'b--','lineWidth',lw);

legend([ab,ar,ag],{'ana','bru','trm'});
xlim([0 1.5])
ylim([0 max(pdf_cr_trm)*1.1])

