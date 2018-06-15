clear;
addpath(fullfile(libgit,'ext-ddm'));
% projectdir = '/hdd/Cloud/Research/matlab/exp/dm-opensimon/v06';
projectdir = '/hdd/Cloud/Research/matlab/exp/dm-sz2/analysis/v19';
f = 'sz_eeg';
cd(projectdir);
%%
dt = fullfile('sim',f);
d = dir(dt);d = {d.name};d = d(contains(d,f));
for ix_d = 1:length(d)
    clear sr;
    fprintf('%s\n',d{ix_d});
    sr = load(fullfile(dt,d{ix_d}));
    sr = sr.obj;
    
    sr.opt.dodoublefit = false;
        
    
%     if any(contains(sr.ddm_print_search,'z_'))
%         sr.ddm_delete_saved('-f');
%     end
    sr.s.reinit = false;
    
    sr.ddm_fit_init(sr.fit(end).p);
    sr.ddm_fit;
%     
    improvement = sr.fit(end-1).nll - sr.fit(end).nll;
    forcesave = false;
    
    th = 2e-3;
    if (improvement>th)
        fprintf('%s improved!!! Saving...\n',ix_d);
    end
    if (improvement>th)||forcesave
        sr.ddm_save;
    end
end