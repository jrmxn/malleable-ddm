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
        
    sr.ddm_fit_init(sr.fit(end).p);
    
    sr.s.reinit = false;
    sr.ddm_fit;
    
    improvement = sr.fit(end-1).nll - sr.fit(end).nll;
    if improvement>2e-3
        fprintf('%s improved!!! Saving...\n',ix_d);
        sr.ddm_save;
    end
end