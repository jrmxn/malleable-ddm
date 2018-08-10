function test_ddm_run_conflict(model_type,varargin)
%% Example function for fitting the SE-DDM
% example use:
% test_ddm_run_conflict('SE-DDM')
%%
d.vec_sub = 1:3;
v = inputParser;
addOptional(v,'vec_sub',d.vec_sub);
parse(v,varargin{:})
v = v.Results;
d = [];clear d;
%%
addpath(fullfile('..','auxf'));
addpath(fullfile('..'));
path_data = fullfile('testing.csv');
n_fits = 5;
vec_sub = v.vec_sub;
%% Fit a basic DDM to get in the ballpark
sr_loaded_fit = cell(1,length(vec_sub));
for ix_sub = 1:length(vec_sub)
    clear sr;
    rng(vec_sub(ix_sub));
    %Set sr to be a ddm_def_conflict object:
    sr = ddm_def_conflict;
    sr.subject = sprintf('sub%02d',vec_sub(ix_sub));
    sr.path_data = path_data;
    sr.info.description = 'seed';
    sr.extra_string = '_seed';
    % get the key to index the model parameters
    mk = sr.ddm_get_instance('keyr');
    % set the search parameters for the model
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    %set the model parameters (the same as the search parameters plus two
    %more)
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    
    sr.ddm_init(id_model_de,id_search_de);
    
    if not(exist(sr.ddm_get_save_path,'file')==2)
        sr.ddm_pdf = @(a,b) sr.ddm_prt_ana(a,b);
        sr.s.reinit = true;
        sr.ddm_fit;
        sr.ddm_save;
    else
        sr = load(sr.ddm_get_save_path);
        sr = sr.obj;
    end
    sr_loaded_fit{1,vec_sub(ix_sub)} = sr.fit(end);
end
%% full model
fprintf('Fitting full model...\n');
for ix_sub = 1:length(vec_sub)
    clear sr;
    rng(vec_sub(ix_sub));r = randi(1e4);
    sr = ddm_def_conflict;
    sr.subject = sprintf('sub%02d',vec_sub(ix_sub));
    disp(sr.subject);
    sr.path_data = path_data;
    sr.info.description = model_type;
    mk = sr.ddm_get_instance('keyr');
    %
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    
    if strcmpi(model_type,'SE-DDM')
        id_search(mk.zc) = 1;
        id_search(mk.b) = 1;
    elseif strcmpi(model_type,'attention')
        id_search(mk.b) = 1;
    elseif strcmpi(model_type,'bias')
        id_search(mk.zc) = 1;
    elseif strcmpi(model_type,'DDM')
    else
        error('Wrong model type');
    end
    
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    
    sr.ddm_init(id_model_de,id_search_de);
    
    %now that we are actually fitting models that potentially have non
    %standard DDM terms, change to the transition matrix approach:
    sr.ddm_pdf = @(a,b) sr.ddm_pdf_trm(a,b,sr.s.dx);
    
    %this is just for resuming
    if (exist(sr.ddm_get_save_path,'file')==2)
        sr_load = load(sr.ddm_get_save_path);
        sr_load = sr_load.obj;
        n_done = length(sr_load.fit);
        [~,best_ix] = min(cellfun(@(x) x.nll,{sr_load.fit}));
        best_fit = sr_load.fit(best_ix);
        clear sr;
        sr = sr_load;
        clear sr_load;
    else
        n_done = 0;
    end
    
    for ix_fit = n_done+1:n_fits
        fprintf('Fit index: %d\n',ix_fit)
        rng(vec_sub(ix_sub) + r*(ix_fit-1));
        if ix_fit == 1
            % for the first fit use the optimum we had found with the
            % simple fit
            sr.s.reinit = false;
            sr.ddm_fit_init(sr_loaded_fit{1,vec_sub(ix_sub)}.p);
        elseif (ix_fit == n_fits)
            % for the last fit, take the best so far and repeat at high
            % resolution:
            sr.s.reinit = false;
            sr.ddm_fit_init(best_fit.p);
            sr.s.dx = 0.005;
            sr.opt.dodoublefit = false;
        else
            % in general, just redraw the initialisation parameters
            sr.s.reinit = true;
        end
        sr.ddm_fit;
        sr.ddm_save;
    end
end
end