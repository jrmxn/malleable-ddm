clear;


f_path_data = 'testing.csv';
n_subs = 21;
f_savepath = cell(1,n_subs);
for ix_sub = 19:n_subs
    clear sr;
    sr = ddm_def_conflict;
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
    
%     sr.subject = subject;
%     sr.path_data = f_path_data;
    
    sr.ddm_init(id_model,id_search);
    sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
    sr.ddm_fit;
    sr.ddm_fit;
    sr.ddm_save;
    px = sr.fit(end).p;px.c = 1;
    %make 100 draws from the pdf
%     sr.ddm_data_draw(px,100);
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
    f_savepath{ix_sub} = sr.ddm_save;
end
for ix_sub = 1:n_subs
    clearvars -except f_savepath;
    sr = load(f_savepath);sr = sr.obj;
    tic;
    sr.ddm_mcmc('mccount',100e3);
    toc;
    sr.ddm_save;
end