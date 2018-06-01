clear;
addpath('..');

niter = 1000;
% f_path_data = 'testing.csv';

for ix_sub = 1:21
    clear sr;
    sr = ddm_def_conflict_h;
    sr.subject = sprintf('sub%02d',ix_sub);
    mk = sr.ddm_get_instance('keyr');
    %%
    id_model = sr.debi_model(0,'de','bi');
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model(mk.a) = 1;
    id_model(mk.t) = 1;
    id_model(mk.v) = 1;
    id_model(mk.st) = 1;
    id_model = sr.debi_model(id_model,'bi','de');
    
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search = sr.debi_model(id_search,'bi','de');
    
    
    sr.ddm_init(id_model,id_search);
    sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
    sr.ddm_fit;
    sr.ddm_fit;
    f_savepath = sr.ddm_save;
    %% Test to use previous model as start point
    % clearvars -except f_savepath;
    % sr = load(f_savepath);sr = sr.obj;
    mk = sr.ddm_get_instance('keyr');
    id_model = sr.debi_model(sr.id_model,'de','bi');
    id_model(mk.zc) = 1;
    id_model(mk.b) = 1;
    id_model = sr.debi_model(id_model,'bi','de');
    id_search = sr.debi_model(sr.id_search,'de','bi');
    id_search(mk.zc) = 1;
    id_search(mk.b) = 1;
    id_search = sr.debi_model(id_search,'bi','de');
    %re-initiate the model
    sr.ddm_init(id_model,id_search);
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    sr.ddm_fit;
    sr.ddm_fit;
    % f_savepath_base_conflict = sr.ddm_save;
    %%
    % f_savepath_base_conflict = 'sim/conflict_h/sub01_x4112384x_x2015232x_conflict_h.mat';
    %% Test to use previous model as start point
    sr_bh = copy(sr);
    mk = sr_bh.ddm_get_instance('keyr');
    id_model = sr_bh.debi_model(sr_bh.id_model,'de','bi');
    id_model(mk.bch) = 1;
    id_model(mk.bchn) = 1;
    id_model = sr_bh.debi_model(id_model,'bi','de');
    id_search = sr_bh.debi_model(sr_bh.id_search,'de','bi');
    id_search(mk.bch) = 1;
    id_search(mk.bchn) = 1;
    id_search = sr_bh.debi_model(id_search,'bi','de');
    %re-initiate the model
    sr_bh.ddm_init(id_model,id_search);
    sr_bh.opt.MaxIter = niter;
    sr_bh.ddm_fit;
    sr_bh.ddm_fit;
    sr_bh.ddm_save;
    clear sr_bh;
    %%
    % clearvars -except f_savepath_base_conflict;
    % sr = load(f_savepath_base_conflict);sr = sr.obj;
    sr_vh = copy(sr);
    mk = sr_vh.ddm_get_instance('keyr');
    id_model = sr_vh.debi_model(sr_vh.id_model,'de','bi');
    id_model(mk.vch) = 1;
    id_model(mk.vchn) = 1;
    id_model = sr_vh.debi_model(id_model,'bi','de');
    id_search = sr_vh.debi_model(sr_vh.id_search,'de','bi');
    id_search(mk.vch) = 1;
    id_search(mk.vchn) = 1;
    id_search = sr_vh.debi_model(id_search,'bi','de');
    %re-initiate the model
    sr_vh.ddm_init(id_model,id_search);
    sr_vh.opt.MaxIter = niter;
    sr_vh.ddm_fit;
    sr_vh.ddm_fit;
    sr_vh.ddm_save;
    clear sr_vh;
    %%m
    % clearvars -except f_savepath_base_conflict;
    % sr = load(f_savepath_base_conflict);sr = sr.obj;
    sr_xbh = copy(sr);
    mk = sr_xbh.ddm_get_instance('keyr');
    id_model = sr_xbh.debi_model(sr_xbh.id_model,'de','bi');
    id_model(mk.zch) = 1;
    id_model(mk.zchn) = 1;
    id_model = sr_xbh.debi_model(id_model,'bi','de');
    id_search = sr_xbh.debi_model(sr_xbh.id_search,'de','bi');
    id_search(mk.zch) = 1;
    id_search(mk.zchn) = 1;
    id_search = sr_xbh.debi_model(id_search,'bi','de');
    %re-initiate the model
    sr_xbh.ddm_init(id_model,id_search);
    sr_xbh.opt.MaxIter = niter;
    sr_xbh.ddm_fit;
    sr_xbh.ddm_fit;
    sr_xbh.ddm_save;
    clear sr_xbh;
end