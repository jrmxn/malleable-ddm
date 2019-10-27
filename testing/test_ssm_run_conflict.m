function test_ssm_run_conflict(model_type,varargin)
%% Example function for fitting the SE-ssm
% example use:
% test_ssm_run_conflict('SE-ssm')
%%
d.vec_sub = 1:3;
v = inputParser;
addOptional(v,'vec_sub',d.vec_sub);
parse(v,varargin{:})
v = v.Results;
d = [];clear d;
%%
addpath(fullfile('..'));
addpath(fullfile('..','ext_mod'));
path_data = fullfile('testing.csv');
n_fits = 5;
vec_sub = v.vec_sub;
%% Fit a basic ssm to get in the ballpark
sr_loaded_fit = cell(1,length(vec_sub));
for ix_sub = 1:length(vec_sub)
    clear sr;
    rng(vec_sub(ix_sub));
    %Set sr to be a ssm_def_conflict object:
    sr = ssm_def_conflict;
    sr.subject = sprintf('sub%02d',vec_sub(ix_sub));
    sr.path_data = path_data;
    sr.info.description = 'seed';
    sr.extra_string = '_seed';
    % get the key to index the model parameters
    mk = sr.ssm_get_instance('keyr');
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
    
    sr.ssm_init(id_model_de,id_search_de);
    
    if not(exist(sr.ssm_get_save_path,'file')==2)
        sr.ssm_pdf = @(a,b) sr.ssm_prt_ana(a,b);
        sr.s.reinit = true;
        sr.ssm_fit;
        sr.ssm_save;
    else
        sr = load(sr.ssm_get_save_path);
        sr = sr.obj;
    end
    sr_loaded_fit{1,vec_sub(ix_sub)} = sr.fit(end);
end
%% full model
fprintf('Fitting full model...\n');
for ix_sub = 1:length(vec_sub)
    clear sr;
    rng(vec_sub(ix_sub));r = randi(1e4);
    sr = ssm_def_conflict;
    sr.subject = sprintf('sub%02d',vec_sub(ix_sub));
    disp(sr.subject);
    sr.path_data = path_data;
    sr.info.description = model_type;
    mk = sr.ssm_get_instance('keyr');
    %
    id_search = sr.debi_model(0,'de','bi');
    id_search(mk.a) = 1;
    id_search(mk.t) = 1;
    id_search(mk.v) = 1;
    id_search(mk.st) = 1;
    
    if strcmpi(model_type,'SE-ssm')
        id_search(mk.zc) = 1;
        id_search(mk.b) = 1;
    elseif strcmpi(model_type,'attention')
        id_search(mk.b) = 1;
    elseif strcmpi(model_type,'bias')
        id_search(mk.zc) = 1;
    elseif strcmpi(model_type,'ssm')
    else
        error('Wrong model type');
    end
    
    id_search_de = sr.debi_model(id_search,'bi','de');
    
    id_model = id_search;
    id_model(mk.s) = 1;
    id_model(mk.z) = 1;
    id_model_de = sr.debi_model(id_model,'bi','de');
    
    sr.ssm_init(id_model_de,id_search_de);
    
    %now that we are actually fitting models that potentially have non
    %standard ssm terms, change to the transition matrix approach:
    sr.ssm_pdf = @(a,b) sr.ssm_pdf_trm(a,b,sr.s.dx);
    
    %this is just for resuming
    if (exist(sr.ssm_get_save_path,'file')==2)
        sr_load = load(sr.ssm_get_save_path);
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
            sr.ssm_fit_init(sr_loaded_fit{1,vec_sub(ix_sub)}.p);
        elseif (ix_fit == n_fits)
            % for the last fit, take the best so far and repeat at high
            % resolution:
            sr.s.reinit = false;
            sr.ssm_fit_init(best_fit.p);
            sr.s.dx = 0.005;
            sr.opt.dodoublefit = false;
        else
            % in general, just redraw the initialisation parameters
            sr.s.reinit = true;
        end
        sr.ssm_fit;
        sr.ssm_save;
    end
end
end