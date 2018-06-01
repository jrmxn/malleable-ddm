clear;
addpath('..');
N = 8;
X = nan(N,8,3);
for ix_sub = 1:N
    f = 'base_ana';
    % f = 'base_ana_prt';
    % f = 'base_trm';sub01_x4079616x_x1966080x_base_trm
    obj = load(sprintf('/hdd/Local/gitprojects/ext-ddm/sim/%s/sub%02d_x4079616x_x1966080x_%s.mat',f,ix_sub,f));
    sr = obj.obj;
    X(ix_sub,:,1) = struct2array(sr.fit(end).p);
    %%
    % f = 'base_ana';
    f = 'base_ana_prt';
    % f = 'base_trm';
    obj = load(sprintf('/hdd/Local/gitprojects/ext-ddm/sim/%s/sub%02d_x4079616x_x1966080x_%s.mat',f,ix_sub,f));
    sr = obj.obj;
    X(ix_sub,:,2) = struct2array(sr.fit(end).p);
    %%
    % f = 'base_ana';
    % f = 'base_ana_prt';
    f = 'base_trm';
    obj = load(sprintf('/hdd/Local/gitprojects/ext-ddm/sim/%s/sub%02d_x4079616x_x1966080x_%s.mat',f,ix_sub,f));
    sr = obj.obj;
    X(ix_sub,:,3) = struct2array(sr.fit(end).p);
end