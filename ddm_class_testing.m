clear;
sr = ddm_def('');
mk = sr.ddm_def_instance('keyr');
%%
id_model = sr.debi_model(0,'de','bi');
id_model(mk.s) = 1;
id_model(mk.a) = 1;
id_model(mk.t) = 1;
id_model(mk.v) = 1;
id_model = sr.debi_model(id_model,'bi','de');

id_fit = sr.debi_model(0,'de','bi');
id_fit(mk.s) = 0;
id_fit(mk.a) = 1;
id_fit(mk.t) = 1;
id_fit(mk.v) = 1;
id_fit = sr.debi_model(id_fit,'bi','de');

sr.subject = 'sub01';
sr.path_data = 'testing.csv';

sr = sr.ddm_init(id_model,id_fit);
%If this isn't run manually it runs automatically.
%But it allows modifications to be made for new branches of model fits - 
%E.g. initialise a new model, but then get the initial p from somewhere
%else. Or write other search initialisation functions.
sr = sr.ddm_fit;
sr = sr.ddm_fit_init;
sr = sr.ddm_fit;
sr.ddm_save;
%%
px = sr.fit(end).p;px.c = 1;
%make 100 draws from the pdf
sr.ddm_data_draw(px,100);
%% Test to use previous model as start point
clear;
sr = load('sub01_x3932160x_x1835008x.mat');sr = sr.obj;
mk = sr.ddm_def_instance('keyr');
id_model = sr.debi_model(sr.id_model,'de','bi');
id_model(mk.xb) = 1;
id_model(mk.b) = 1;
id_model = sr.debi_model(id_model,'bi','de');
id_fit = sr.debi_model(sr.id_fit,'de','bi');
id_fit(mk.xb) = 1;
id_fit(mk.b) = 1;
id_fit = sr.debi_model(id_fit,'bi','de');
%re-initiate the model
sr = sr.ddm_init(id_model,id_fit);
sr = sr.ddm_fit;
sr.ddm_save;
%%
clear;
sr = load('sub01_x4128768x_x2031616x.mat');sr = sr.obj;
sr = sr.ddm_mcmc('mccount',100);

%%
% figure(2);
% [C,lags,ESS]=eacorr(sr.mcmc.models);
% plot(lags,C,'.-',lags([1 end]),[0 0],'k');
% grid on
% xlabel('lags')
% ylabel('autocorrelation');
% text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
% title('Markov Chain Auto Correlation')