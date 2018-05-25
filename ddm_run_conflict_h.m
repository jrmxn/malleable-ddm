clear;
sr = ddm_def_conflict_h;
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
/home/mcintosh
sr.subject = 'sub01';
sr.path_data = 'testing.csv';

sr.ddm_init(id_model,id_search);
% sr.opt.MaxIter = 10;
sr.ddm_fit;
sr.ddm_fit;
f_savepath = sr.ddm_save;
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
f_savepath_base_conflict = sr.ddm_save;
%%
f_savepath_base_conflict = 'sim/conflict_h/sub01_x4112384x_x2015232x_conflict_h.mat';
%% Test to use previous model as start point
clearvars -except f_savepath_base_conflict;
sr = load(f_savepath_base_conflict);sr = sr.obj;
mk = sr.ddm_get_instance('keyr');
id_model = sr.debi_model(sr.id_model,'de','bi');
id_model(mk.bhc) = 1;
id_model(mk.bhnc) = 1;
id_model = sr.debi_model(id_model,'bi','de');
id_search = sr.debi_model(sr.id_search,'de','bi');
id_model(mk.bhc) = 1;
id_model(mk.bhnc) = 1;
id_search = sr.debi_model(id_search,'bi','de');
%re-initiate the model
sr.ddm_init(id_model,id_search);
% sr.opt.MaxIter = 10;
sr.ddm_fit;
sr.ddm_fit;
%%
clearvars -except f_savepath_base_conflict;
sr = load(f_savepath_base_conflict);sr = sr.obj;
mk = sr.ddm_get_instance('keyr');
id_model = sr.debi_model(sr.id_model,'de','bi');
id_model(mk.vhc) = 1;
id_model(mk.vhnc) = 1;
id_model = sr.debi_model(id_model,'bi','de');
id_search = sr.debi_model(sr.id_search,'de','bi');
id_model(mk.bhc) = 1;
id_model(mk.bhnc) = 1;
id_search = sr.debi_model(id_search,'bi','de');
%re-initiate the model
sr.ddm_init(id_model,id_search);
% sr.opt.MaxIter = 10;
sr.ddm_fit;
sr.ddm_fit;
%%
clearvars -except f_savepath_base_conflict;
sr = load(f_savepath_base_conflict);sr = sr.obj;
mk = sr.ddm_get_instance('keyr');
id_model = sr.debi_model(sr.id_model,'de','bi');
id_model(mk.xbhc) = 1;
id_model(mk.xbhnc) = 1;
id_model = sr.debi_model(id_model,'bi','de');
id_search = sr.debi_model(sr.id_search,'de','bi');
id_model(mk.bhc) = 1;
id_model(mk.bhnc) = 1;
id_search = sr.debi_model(id_search,'bi','de');
%re-initiate the model
sr.ddm_init(id_model,id_search);
% sr.opt.MaxIter = 10;
sr.ddm_fit;
sr.ddm_fit;