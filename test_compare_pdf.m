clear all;
%%
sr = load('sim/conflict/sub01_x4141056x_x0995328x_conflict.mat');sr = sr.obj;
p = sr.fit(end).p;
sr.s.dt = 1e-3;
sr.s.nits = 100e3;
sr.s.dx = 5e-3;

lt = [0:1e-3:2];
p.c = 1;
% p.a = p.a;
% p.v = 2;
% p.st = 0.0;
% p.sv = 0.0;
% p.sz = 0.0;
% p.z = 0.5;
tic
[pdf_cw_bru,pdf_cr_bru,rt_bru,cdf_dow,cdf_ups] = sr.ddm_pdf_bru(p,lt,sr.s.nits);
toc
%%
p2 = p;
%%
tic
[pdf_cw_ana,pdf_cr_ana,rt_ana] = sr.ddm_pdf_ana(p2,lt);
toc
%%
tic
[pdf_cw_trm,pdf_cr_trm,rt_trm] = sr.ddm_pdf_trm(p,lt,sr.s.dx);
toc
%%
close all;
clf;
hold on;
% subplot(2,1,1);cla;hold on;
% subplot(2,1,2);cla;hold on;
% full_pdf(0.3,p.v,0,p.a,0 ,0 ,p.t,p.st,err)
lw = 1;

ab = plot(rt_ana,pdf_cw_ana,'b','lineWidth',lw);
plot(rt_ana,pdf_cr_ana,'b--','lineWidth',lw);

ar = plot(rt_bru,pdf_cr_bru,'r','lineWidth',lw);
plot(rt_bru,pdf_cw_bru,'r','lineWidth',lw);

ag = plot(rt_trm,pdf_cw_trm,'g','lineWidth',lw);
plot(rt_trm,pdf_cr_trm,'g','lineWidth',lw);

legend([ab,ar,ag],{'ana','bru','trm'});
xlim([0 1.5])
ylim([0 max(pdf_cr_ana)*1.1])

% p3 = p2;p3.sz = 0;
% [pdf_cw_ana,pdf_cr_ana,rt_ana] = sr.ddm_pdf_ana(p3,lt);
% plot(rt_ana,pdf_cr_ana,'k--','lineWidth',lw);