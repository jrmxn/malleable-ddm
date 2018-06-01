clear;
f_path_data = 'testing.csv';
addpath('..');
for ix_fit_type = [1,2,3]
    for ix_sub = 1:3
        
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
        id_model(mk.sz) = 1;
%         id_model(mk.sv) = 1;
        id_model_de = sr.debi_model(id_model,'bi','de');
        
        id_search = sr.debi_model(0,'de','bi');
        id_search(mk.s) = 0;
        id_search(mk.z) = 0;
        id_search(mk.a) = 1;
        id_search(mk.t) = 1;
        id_search(mk.v) = 1;
        id_search(mk.st) = 1;
        id_search(mk.sz) = 1;
%         id_search(mk.sv) = 1;
        
        id_search_de = sr.debi_model(id_search,'bi','de');

        sr.ddm_init(id_model_de,id_search_de);
        if ix_fit_type == 2
            sr.modelclass = 'base_ana';
            sr.ddm_pdf = @(a,b) sr.ddm_pdf_ana(a,b);
        elseif ix_fit_type ==1
            sr.modelclass = 'base_ana_prt';
            sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
        elseif ix_fit_type==3
            sr.modelclass = 'base_trm';
            sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
        end
        sr.ddm_fit;
        sr.ddm_fit;
        sr.ddm_save;
    end
end