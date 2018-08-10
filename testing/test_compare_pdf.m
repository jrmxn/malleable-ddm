clear;
addpath('..');
%%
sr = load('sim/base_ana_prt/sub01_x543313362944x_x130996502528x_base_ana_prt.mat');sr = sr.obj;
p = sr.fit(end).p;
sr.s.dt = 1e-3;
sr.s.nits = 100e3;
sr.s.dx = 5e-3;

lt = [0:1e-3:5];
p.t = 0.25;
p.s = 1;
p.c = 0;
p.v = 1.5;
p.st = 0.0;
p.sv = 0.0;
p.sz = 0.0;
p.z = 0.5;
p.a = 1.0;

tic
% p.v = 1.5;
[pdf_cw_bru,pdf_cr_bru,rt_bru,cdf_dow,cdf_ups] = sr.ddm_pdf_bru(p,lt,sr.s.nits);
toc
%%
tic
% p.v = 1.5;
[pdf_cw_trm,pdf_cr_trm,rt_trm] = sr.ddm_pdf_trm(p,lt,sr.s.dx);
toc
%%
close all;clf;hold on;
lw = 1;

ar = plot(rt_bru,pdf_cr_bru,'r','lineWidth',lw);
plot(rt_bru,pdf_cw_bru,'--r','lineWidth',lw);

ag = plot(rt_trm,pdf_cr_trm,'g','lineWidth',lw);
plot(rt_trm,pdf_cw_trm,'--g','lineWidth',lw);

legend([ar,ag],{'bru','trm'});
xlim([0 1.5])
ylim([0 max(pdf_cr_trm)*1.1])