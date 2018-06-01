clear;
addpath('..');

t = 3;
N_sub = 20;
X = nan(N_sub,2);
for ix_sub = 1:N_sub
    if t == 1
        f = sprintf('sub%02d_x4113152x_x2016000x_conflict_h',ix_sub);
    elseif t == 2
        f = sprintf('sub%02d_x4115456x_x2018304x_conflict_h',ix_sub);
    elseif t ==3
        f = sprintf('sub%02d_x4124672x_x2027520x_conflict_h',ix_sub);
    end
    
    sr = load(sprintf('/hdd/Local/gitprojects/dm-conflictmodel/sim/conflict_h/%s.mat',f));
    sr = sr.obj;
    if t==1
        X(ix_sub,:) = [sr.fit(end).p.xbhc,sr.fit(end).p.xbhnc];
    elseif t==2
        X(ix_sub,:) = [sr.fit(end).p.vhc,sr.fit(end).p.vhnc];
    elseif t==3
        X(ix_sub,:) = [sr.fit(end).p.bhc,sr.fit(end).p.bhnc];
    end
end

if t==1
    fprintf('xbhc\txbhnc\n');
elseif t==2
    fprintf('vhc\tvhnc\n');
elseif t==3
    fprintf('bhc\tbhnc\n');
end
X