clear;
addpath('..');
[sub_data_struc,sub_cell,sz] = base_def_subjects('subjects','ALL');
path_data = fullfile(pwd,'testing_sz.csv');
%%
for ix_sub = 1:length(sub_cell)
    
    clear sr;
    sr = ddm_def_sz;
    sr.subject = sub_cell{ix_sub};
    sr.path_data = path_data;
    %     sr.get_data;% in this case our instance is dependent on data so get it first.
    mk = sr.ddm_get_instance('keyr');
    %% model def
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.vm3) = 1;id_search(mk.vm4) = 1;id_search(mk.vm5) = 1;
    id_search(mk.vz) = 1;
    id_search(mk.vp3) = 1;id_search(mk.vp4) = 1;id_search(mk.vp5) = 1;
    if contains(sr.subject,sz)
        id_search(mk.vm1) = 1;id_search(mk.vm2) = 1;
        id_search(mk.vp1) = 1;id_search(mk.vp2) = 1;
    end
    id_search(mk.st) = 1;
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    %%
    sr.ddm_init(id_model_de,id_search_de);
    sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
    
    sr.ddm_fit;
    sr.ddm_fit;
    sr.ddm_save;
end