clear;
sub_cell{1} = 'AHLAJ';sz(1) = true;
sub_cell{2} = 'QQYAK';sz(2) = false;
for ix_sub = 1:length(sub_cell)
    
    clear sr;
    sr = ddm_def_sz_base;
    sr.subject = sub_cell{ix_sub};
%     sr.get_data;% in this case our instance is dependent on data so get it first.
    mk = sr.ddm_get_instance('keyr');
    %%
    id_model = sr.debi_model(0,'de','bi');
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model(mk.a) = 1;
    id_model(mk.t) = 1;
    id_model(mk.vm1) = 1;id_model(mk.vm2) = 1;id_model(mk.vm3) = 1;id_model(mk.vm4) = 1;id_model(mk.vm5) = 1;
    id_model(mk.vz) = 1;
    id_model(mk.vp1) = 1;id_model(mk.vp2) = 1;id_model(mk.vp3) = 1;id_model(mk.vp4) = 1;id_model(mk.vp5) = 1;
    id_model(mk.st) = 1;
    id_model = sr.debi_model(id_model,'bi','de');
    
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.s) = 0;
    id_search(mk.z) = 0;
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.vm3) = 1;id_search(mk.vm4) = 1;id_search(mk.vm5) = 1;
    id_search(mk.vz) = 1;
    id_search(mk.vp3) = 1;id_search(mk.vp4) = 1;id_search(mk.vp5) = 1;
    if not(sz)
        id_search(mk.vm1) = 1;id_search(mk.vm2) = 1;
        id_search(mk.vp1) = 1;id_search(mk.vp2) = 1;
    end
    id_search(mk.st) = 1;
    
    id_search = sr.debi_model(id_search,'bi','de');
    
    sr.ddm_init(id_model,id_search);
    
    sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
    
    sr.ddm_fit;
    sr.ddm_fit;
    sr.ddm_save;
end