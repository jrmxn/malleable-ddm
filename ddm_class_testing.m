clear;
sr = ddm_def('');
mk = sr.get_modeldef('keyr');
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
sr.opt.MaxIter = 200;
%If this isn't run manually it runs automatically.
%But it allows modifications to be made for new branches of model fits - 
%E.g. initialise a new model, but then get the initial p from somewhere
%else. Or write other search initialisation functions.
sr = sr.ddm_search;
sr = sr.ddm_search_init;
sr = sr.ddm_search;
sr.ddm_save;

px = sr.fit(end).p;px.c = 1;
%make 100 draws from the pdf
sr.ddm_data_draw(px,100);
%%
clear;
sr = load('sub01_x3932160x_x1835008x.mat');sr = sr.obj;