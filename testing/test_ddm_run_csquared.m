% function test_ddm_run_conflict
n_subs = 21;%%
vec_sub = 1:n_subs
% :n_subs;
% clearvars -except vec_sub vec_sub_oi;
% sub_data_struc = base_def_subjects('subjects','ALL');
% which_ddm_def = which('ddm_def');
% if isempty(which_ddm_def)
addpath('..');
% addpath(fullfile(libgit,'dm_ext_ddm'));
% end
path_data = fullfile(pwd,'testing.csv');
% path_data = fullfile(pwd,'local_data','20180531_simon_beh.csv');


n_subs = 21;
%% Fit a basic DDM to get in the ballpark
for ix_sub = 1:length(vec_sub)
    clear sr;
    rng(ix_sub);
    sr = ddm_def_csquared;
    sr.subject = sprintf('sub%02d',ix_sub);
    sr.path_data = path_data;
    mk = sr.ddm_get_instance('keyr');
    %
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    
    sr.ddm_init(id_model_de,id_search_de);
    
    if not(exist(sr.ddm_get_save_path,'file')==2)
        sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
        sr.s.reinit = true;
        sr.opt.MaxIter = 2500;
        sr.ddm_fit;
        sr.ddm_fit;
%         sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
%         sr.s.reinit = false;
%         sr.ddm_fit;
sr.ddm_save;
    end
end
%% full model
for ix_sub = 1:length(vec_sub)
    clear sr;
    rng(ix_sub);
    sr_load = ddm_def_csquared;
    sr_load.subject = sprintf('sub%02d',ix_sub);
    sr_load.path_data = path_data;
    mk = sr_load.ddm_get_instance('keyr');
    %
    id_search = sr_load.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search_de = sr_load.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr_load.debi_model(id_model,'bi','de');
    
    sr_load.ddm_init(id_model_de,id_search_de);
    
    sr = load(sr_load.ddm_get_save_path);
    sr = sr.obj;clear sr_load;
    
    id_search(mk.zc) = 1;
    id_search(mk.b) = 1;
    id_search_de = sr.debi_model(id_search,'bi','de');
    id_model(mk.zc) = 1;
    id_model(mk.b) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    sr.ddm_init(id_model_de,id_search_de);
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    
    if not(exist(sr.ddm_get_save_path,'file')==2)
        sr.s.reinit = false;
        sr.opt.MaxIter = sr.opt.MaxIter*3;
        sr.ddm_fit;
        sr.ddm_fit;
        sr.ddm_save;
    end
    
end
% for ix_sub = 1:length(vec_sub)
%     clearvars -except f_savepath;
%     sr = load(f_savepath);sr = sr.obj;
%     tic;
%     sr.ddm_mcmc('mccount',100e3);
%     toc;
%     sr.ddm_save;
% end
% rmpath(fullfile(libgit,'dm_ext_ddm'));