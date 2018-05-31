clear;
f_path_data = 'testing.csv';
for ix_sub = 1
    
    clear sr;
    sr = ddm_def_base;
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
    id_search(mk.s) = 0;
    id_search(mk.z) = 0;
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    
    id_search = sr.debi_model(id_search,'bi','de');
    
    sr.ddm_init(id_model,id_search);
        sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
    sr.ddm_fit;
%    
    sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
    [nll_prt_ana,~,~,~,ll_vec_prt_ana] = ddm_cost_pdf_nll(sr,[],sr.fit.p);
    
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_ana(a,b);
    [nll_pdf_ana,~,~,~,ll_vec_pdf_ana] = ddm_cost_pdf_nll(sr,[],sr.fit.p);
    
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    [nll_pdf_trm,~,~,~,ll_vec_pdf_trm] = ddm_cost_pdf_nll(sr,[],sr.fit.p);

    sr.ddm_pdf = @(a,b) sr.ddm_pdf_bru(a,b,25e3);
    [nll_pdf_bru,~,~,~,ll_vec_pdf_bru] = ddm_cost_pdf_nll(sr,[],sr.fit.p);
    %
    [nll_prt_ana,nll_pdf_ana,nll_pdf_trm,nll_pdf_bru]
    T = [ll_vec_prt_ana,ll_vec_pdf_ana,ll_vec_pdf_trm,ll_vec_pdf_bru];
    T = array2table(T);
    T.Properties.VariableNames = {'prt_ana','pdf_ana','pdf_trm','pdf_bru'};
    Tl = T.Properties.VariableNames;
    for ix_1 = 1:length(Tl)
        for ix_2 = 1:length(Tl)
            subplot(length(Tl),length(Tl), length(Tl)*(ix_2-1) + ix_1 )
            plot(T.(Tl{ix_1}),T.(Tl{ix_2}),'.')
            xlabel(strrep(Tl{ix_1},'_',' '))
            ylabel(strrep(Tl{ix_2},'_',' '))
        end
    end
end
