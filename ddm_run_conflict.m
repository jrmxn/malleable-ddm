clear;
sr = ddm_def_conflict;
mk = sr.ddm_get_instance('keyr');
%%
id_model = sr.debi_model(0,'de','bi');
id_model(mk.s) = 1;
id_model(mk.a) = 1;
id_model(mk.t) = 1;
id_model(mk.v) = 1;
id_model(mk.st) = 1;
id_model = sr.debi_model(id_model,'bi','de');

id_search = sr.debi_model(0,'de','bi');
id_search(mk.s) = 0;
id_search(mk.a) = 1;
id_search(mk.t) = 1;
id_search(mk.v) = 1;
id_search(mk.st) = 1;
id_search = sr.debi_model(id_search,'bi','de');

sr.subject = 'sub01';
sr.path_data = 'testing.csv';

sr.ddm_init(id_model,id_search);
% sr.opt.MaxIter = 10;
sr.ddm_fit;
sr.ddm_fit;
f_savepath = sr.ddm_save;
px = sr.fit(end).p;px.c = 1;
%make 100 draws from the pdf
sr.ddm_data_draw(px,100);
%% Test to use previous model as start point
clearvars -except f_savepath;
sr = load(f_savepath);sr = sr.obj;
mk = sr.ddm_get_instance('keyr');
id_model = sr.debi_model(sr.id_model,'de','bi');
id_model(mk.xb) = 1;
id_model(mk.b) = 1;
id_model = sr.debi_model(id_model,'bi','de');
id_search = sr.debi_model(sr.id_search,'de','bi');
id_search(mk.xb) = 1;
id_search(mk.b) = 1;
id_search = sr.debi_model(id_search,'bi','de');
%re-initiate the model
sr.ddm_init(id_model,id_search);
% sr.opt.MaxIter = 10;
sr.ddm_fit;
sr.ddm_fit;
f_savepath = sr.ddm_save;
%%
clearvars -except f_savepath;
sr = load(f_savepath);sr = sr.obj;
tic;
sr.ddm_mcmc('mccount',50e3);
toc;
sr.ddm_save;
%%
figure(2);
[C,lags,ESS]=eacorr(sr.mcmc(end).models);
plot(lags,C,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
title('Markov Chain Auto Correlation')