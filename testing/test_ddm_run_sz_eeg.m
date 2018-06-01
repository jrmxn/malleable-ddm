clear;
addpath('..');
sub_cell{1} = 'AHLAJ';sz(1) = true;
% sub_cell{2} = 'QQYAK';sz(2) = false;
%%

%%
for ix_sub = 1:length(sub_cell)
    
    clear sr;
    sr = ddm_def_sz_eeg;
    sr.subject = sub_cell{ix_sub};
    fprintf('%s\n',sr.subject);
    %     sr.get_data;% in this case our instance is dependent on data so get it first.
    mk = sr.ddm_get_instance('keyr');
    %% model definition
    eeg_mod(1).param = 'v';
    eeg_mod(1).channel = 'pre_alpha_O';
    
    id_search = sr.debi_model(0,'de','bi');
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
    for ix_eeg_mod = 1:length(eeg_mod)
        ch_str = eeg_mod(ix_eeg_mod).channel;
        p_str = eeg_mod(ix_eeg_mod).param;
        id_search(mk.(sprintf('%s_%s',p_str,ch_str))) = 1;
    end
    id_search_de = sr.debi_model(id_search,'bi','de');
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    %%    
    sr.ddm_init(id_model_de,id_search_de);
        
    sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b,eeg_mod);
    
    sr.ddm_fit;
    sr.ddm_fit;
    sr.ddm_save;
    sr.ddm_mcmc;
    sr.ddm_save;
end