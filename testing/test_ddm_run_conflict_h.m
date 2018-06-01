            
clear;
vec_sub_oi = 1;
n_subs = 21;%%
vec_sub = 1:n_subs;
clearvars -except vec_sub vec_sub_oi;
addpath('..');
% sub_data_struc = base_def_subjects('subjects','ALL');
% path_data = fullfile(pwd,'local_data','20180531_simon_beh.csv');
path_data = fullfile(pwd,'testing.csv');
%%
n_random_start = 1;
sr_loaded_fit = cell(1,length(vec_sub));
for ix_sub = 1:length(vec_sub)
    sr_load = ddm_def_conflict;
    sr_load.subject = sprintf('sub%02d',vec_sub(ix_sub));
    sr_load.path_data = path_data;
    mk = sr_load.ddm_get_instance('keyr');
    %
    id_search = sr_load.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search(mk.zc) = 1;
    id_search(mk.b) = 1;
    id_search_de = sr_load.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr_load.debi_model(id_model,'bi','de');
    
    sr_load.ddm_init(id_model_de,id_search_de);
    
    sr_loaded = load(sr_load.ddm_get_save_path);
    sr_loaded = sr_loaded.obj;
    
    sr_loaded_fit{1,vec_sub(ix_sub)} = sr_loaded.fit(end);
    clear sr_load sr_loaded;
end
%%
fprintf('zch and zchn\n');
for ix_sub = 1:length(vec_sub_oi)
    clear sr;
    sr = ddm_def_conflict_h;
    sr.subject = sprintf('sub%02d',vec_sub_oi(ix_sub));
    sr.path_data = path_data;
    mk = sr.ddm_get_instance('keyr');
    %
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search(mk.zc) = 1;
    id_search(mk.b) = 1;
    %
    id_search(mk.zch) = 1;
    id_search(mk.zchn) = 1;
    %
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    %
    sr.ddm_init(id_model_de,id_search_de);
    sr.ddm_fit_init(sr_loaded_fit{1,vec_sub_oi(ix_sub)}.p);
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    
    if not(exist(sr.ddm_get_save_path,'file')==2)
        sr.s.reinit = false;
        sr.opt.MaxIter = 2500;
        sr.ddm_fit;
        sr.ddm_fit;
    end
    sr.ddm_save;
end
%%
fprintf('bch and bchn\n');
for ix_sub = 1:length(vec_sub_oi)
    clear sr;
    sr = ddm_def_conflict_h;
    sr.subject = sprintf('sub%02d',vec_sub_oi(ix_sub));
    sr.path_data = path_data;
    mk = sr.ddm_get_instance('keyr');
    %
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search(mk.zc) = 1;
    id_search(mk.b) = 1;
    %
    id_search(mk.bch) = 1;
    id_search(mk.bchn) = 1;
    %
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    %
    sr.ddm_init(id_model_de,id_search_de);
    sr.ddm_fit_init(sr_loaded_fit{1,vec_sub_oi(ix_sub)}.p);
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    
    if not(exist(sr.ddm_get_save_path,'file')==2)
        sr.s.reinit = false;
        sr.opt.MaxIter = 2500;
        sr.ddm_fit;
        sr.ddm_fit;
    end
    sr.ddm_save;
end

%%
fprintf('vch and vchn\n');
for ix_sub = 1:length(vec_sub_oi)
    clear sr;
    sr = ddm_def_conflict_h;
    sr.subject = sprintf('sub%02d',vec_sub_oi(ix_sub));
    sr.path_data = path_data;
    mk = sr.ddm_get_instance('keyr');
    %
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search(mk.zc) = 1;
    id_search(mk.b) = 1;
    %
    id_search(mk.vch) = 1;
    id_search(mk.vchn) = 1;
    %
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    %
    sr.ddm_init(id_model_de,id_search_de);
    sr.ddm_fit_init(sr_loaded_fit{1,vec_sub_oi(ix_sub)}.p);
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    
    if not(exist(sr.ddm_get_save_path,'file')==2)
        sr.s.reinit = false;
        sr.opt.MaxIter = 2500;
        sr.ddm_fit;
        sr.ddm_fit;
    end
    sr.ddm_save;
end
rmpath('..')