classdef ddm_def < matlab.mixin.Copyable%instead of handle
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        subject = '';
        modelclass = '';
        id_model = -1+2^4;
        id_search = 1;
        path_data = '';
        fit_ix = 0;
        extra_string = '';
        modelKey = [];
        data = [];
        s = [];
        opt = [];
        fit = [];
        info = [];
        mcmc = [];
        ddm_pdf = [];
        pmod = struct;
    end
    properties (Transient = true)
        objp = [];
    end
    
    methods
        function obj = ddm_def
            %light initialisation so functions can be used easily
            obj.modelclass = '';%used to set file names
            obj.info.version = sprintf('0.0.7');
            obj.info.date = datetime;
            obj.info.description = '';
            try
                %oly writes lowest in the hierarchy, but it's a start
                [ST, I] = dbstack('-completenames', 1);
                obj.info.code = fileread(ST(I).file);
            catch
                warning('Failed to archive code.');
            end
        end
        %         function set.subject(obj,sub)
        % for some reason, set methods don't play nice with minimisation
        % (patternsearch)
        %             fprintf('Setting subject name to %s\n',sub);
        %             obj.subject = sub;
        %         end
        function ddm_writecode(obj,filename)
            fid = fopen(filename,'wt');
            fprintf(fid, '%s',obj.info.code);
            fclose(fid);
        end
        
        function g = ddm_print_search(obj)
            if isempty(obj.modelKey)
                fprintf('Model key not set. Setting it...\n');
                obj.modelKey = obj.ddm_get_instance('keyf');
            end
            g = obj.modelKey(find(obj.debi_model(obj.id_search,'de','bi')));
            disp(g);
        end
        
        
        function ddm_delete_saved(obj,flag)
            if not(exist('flag','var')==1)
                flag = '';
            end
            markfordeletion = false;
            
            if contains(flag,'f')
                markfordeletion = true;
            else
                ip = input('Are you sure you want to delete? [y/n]','s');
                if strcmpi(ip,'y')
                    fprintf('OK. Deleting.\n');
                    markfordeletion = true;
                elseif strcmpi(ip,'n')
                    fprintf('OK. Not deleting.\n');
                else
                    fprintf('What? Not deleteing.\n');
                end
            end
            if markfordeletion
                delete(obj.ddm_get_save_path);
            end
        end
        
        function [p_mat,sr_full] = aux_gather(obj,f_path,id_model_de,id_search_de,sub_cell,minorfin)
            if not(exist('minorfin','var')==1),minorfin = 'fin';end
            if not(isnumeric(id_model_de)&isnumeric(id_search_de))
                error('You put non numeric inputs into aux_gather - expects decimal model specification.');
            end
            obj.id_model = id_model_de;
            obj.id_search = id_search_de;
            obj.subject = '**';
            [~,f_base] = fileparts(obj.ddm_get_save_path);
            f = fullfile(f_path,[f_base,'.mat']);
            vec_empty = [];
            p_mat = struct([]);clear p_mat;%annoying to init properly.
            sr_full = struct([]);clear sr_full;%annoying to init properly.
            
            for ix_sub = 1:length(sub_cell)
                f_ = strrep(f,'**',sub_cell{ix_sub});
                
                if exist(f_,'file') == 2
                    sr = load(f_);
                    sr = sr.obj;
                    sr_full(ix_sub) = sr;
                    if strcmpi(minorfin,'min')
                        [~,ix_min] = min([sr.fit.nll]);
                        sr_fit = sr.fit(ix_min);
                    elseif strcmpi(minorfin,'fin')
                        sr_fit = sr.fit(end);
                    else
                        error('minfin?');
                    end
                    p = sr_fit.p;
                    p.nll = sr_fit.nll;
                    p.aic = sr_fit.aic;
                    p.bic = sr_fit.bic;
                    %p.ll_vec = sr_fit.ll_vec;
                    p.n = sum(not(isnan(sr_fit.ll_vec)));
                    p.subject = sr.subject;
                    p_mat(ix_sub) = p;
                    ix_non_empty = ix_sub;
                else
                    vec_empty = [vec_empty,ix_sub];
                    warning('Missing %s',f_);
                end
            end
            if not(exist('p_mat','var')==1)
                error('No files found, e.g.: %s',f_);
            end
            % fill in the blanks with nans, ''
            % get an example of a subject result that we have
            p_proto = p_mat(ix_non_empty);
            %make sure that p_mat extends to all subjects (even though some
            %are empty)
            p_mat(length(sub_cell)+1) = p_proto;p_mat(length(sub_cell)+1) = [];
            for ix_vec_empty = 1:length(vec_empty)
                p_ = p_mat(vec_empty(ix_vec_empty));
                vec_fn = fieldnames(p_);
                for ix_fn = 1:length(vec_fn)
                    if isnumeric(p_proto.(vec_fn{ix_fn}))
                        v = nan;
                    elseif ischar(p_proto.(vec_fn{ix_fn}))
                        v = '';
                    else
                        error('Undefined case');
                    end
                    p_.(vec_fn{ix_fn}) = v;
                end
                p_mat(vec_empty(ix_vec_empty)) = p_;
            end
            %
            p_mat = struct2table(p_mat);
        end
        
        function ddm_init(obj, id_model,id_search)
            %Set the simulation settings (s), and the optimisation settings
            %(opt). id_model and if_fit (which are decimals,
            % which correspond to binary vectors referencing the
            % parameters) get set.
            % additionally id_search gets translated into xl which holds the
            % fit parameters as strings in a cell.
            if isempty(obj.data)
                get_data(obj);
            end
            obj.modelKey = ddm_get_instance(obj, 'keyf');
            
            obj.id_model = id_model;
            obj.id_search = id_search;
            err_str = 'Bad model spec - check you did not mix up id_model and id_search';
            if not(isnumeric(id_model)&(numel(id_model)==1)),error(err_str);end
            if not(isnumeric(id_search)&(numel(id_search)==1)),error(err_str);end
            
            obj.s.fit_n = sum(ddm_def.debi_model(obj.id_search,'de','bi'));
            obj.s.minAlgo = 'nll';
            obj.s.reinit = false;
            obj.s.dt = 1e-3;
            %             obj.s.ddx = 100;%no longer used
            obj.s.dx = 0.01;%0.01 about = ddx = 100
            %             obj.s.x_bound_scale = 6;
            
            obj.s.T = 5;%this could be set dynamically based on closed form
            obj.s.nits = 50e3;
            %             obj.s.inittype = 'random';
            obj.s.path_data = '';
            obj.ddm_pdf = @(a,b) obj.ddm_pdf_ana(a,b);
            %             obj.ddm_pdf = @(a,b,c) obj.ddm_pdf_bru(a,b,c,obj.s.nits);
            %             obj.ddm_pdf = @(a,b,c) obj.ddm_pdf_trm(a,b,c,obj.s.dx);
            
            id_search_index  = find(ddm_def.debi_model(obj.id_search,'de','bi'));
            
            % Set the parameters over which to optimise from modelType spec
            co = 1;
            for ix_id_search_index = 1:length(id_search_index)
                parameter_string = obj.modelKey{id_search_index(ix_id_search_index)};
                obj.s.xl.(parameter_string) = co;co = co+1;
            end
            
            obj.opt.dodoublefit = true;
            obj.opt.TolX = 1e-5;
            obj.opt.MaxIter = 50e3;
            obj.opt.parallelsearch = true;
            obj.opt.ps_AccelerateMesh = true;%should only do if smooth
            obj.opt.computeAlgo = 'PS';
            if strcmpi(obj.s.minAlgo,'nll')
                obj.opt.h_cost = @obj.ddm_cost_pdf_nll;
            else
                error('minAlgo not defined');
            end
            
        end
        
        function ddm_fit_init(obj,fit_ix_)
            %if you want to run this manually you run it like this:
            %             sr.ddm_fit_init(sr.fit_ix + 1)
            % which then initialises it, but doesn't actually increase sr.fit_ix which
            % is then increased in ddm_fit
            if isempty(obj.data)
                get_data(obj);
            end
            if isstruct(fit_ix_)
                %in the special case where we are inserting a p structure
                %here to initialise...
                init_p_insert = fit_ix_;
                fit_ix_ = 1;%assume
            else
                %in normal cases make sure insert structure is never a
                %match with parameter names
                init_p_insert.x_no_match_x = 0;
            end
            if (fit_ix_==1)||obj.s.reinit
                %if first model fit for this object or we want to reinitialise
                %p is defined from a full set of defaults/random init
                init_p_full_random = obj.ddm_get_instance(['init_' 'random']);
                init_p_full_default = obj.ddm_get_instance(['init_' 'default']);
                % and now we need to reduce it to match the specific model
                % configuration we are testing.
                id_model_index = find(ddm_def.debi_model(obj.id_model,'de','bi'));
                id_search_index  = find(ddm_def.debi_model(obj.id_search,'de','bi'));
                
                for ix_parameter_cell = 1:length(obj.modelKey)
                    parameter_string = obj.modelKey{ix_parameter_cell};
                    if any(strcmp(parameter_string,fieldnames(init_p_insert)))
                        %if we want to insert a p-value structure, then
                        %the inserted p-values get priority
                        init_p_reduced.(parameter_string) = init_p_insert.(parameter_string);
                    elseif any(strcmp(parameter_string,obj.modelKey(id_search_index)))
                        %if we are going to try to optimise this parameter
                        %then get it from a distribution
                        init_p_reduced.(parameter_string) = init_p_full_random.(parameter_string);
                    elseif any(strcmp(parameter_string,obj.modelKey(id_model_index)))
                        %If they are not in search but they are in the
                        %base model def set them to defaults (e.g. usually
                        %                         s = 1, z = 0.5)
                        init_p_reduced.(parameter_string) = init_p_full_default.(parameter_string);
                    else
                        % if the parameter just isn't in this model
                        % definition.
                        init_p_reduced.(parameter_string) = 0;
                    end
                end
                fit_init.p = init_p_reduced;
            else
                fprintf('Resuming from previous fit.\n');
                
                same_model = obj.fit(fit_ix_-1).id_model==obj.id_model;
                same_search = obj.fit(fit_ix_-1).id_search==obj.id_search;
                if not(same_model|same_search),fprintf('With a new model!\n');end
                fit_init.p = obj.fit(fit_ix_-1).p;
            end
            fprintf('Getting initial cost...: ');
            fit_init.nll = obj.opt.h_cost([],fit_init.p);
            fprintf('%0.3f\n',fit_init.nll);
            fit_init.p_lb = obj.ddm_get_instance('lbound');
            fit_init.p_ub = obj.ddm_get_instance('ubound');
            
            obj.fit(fit_ix_).init = fit_init;
        end
        
        function ddm_fit(obj)
            fit_ix_ = obj.fit_ix + 1;
            
            if length(obj.fit)<fit_ix_
                ddm_fit_init(obj,fit_ix_);
            else
                fprintf('Assuming ddm fit was manually initialised.\n');
            end
            
            fit_init = obj.fit(fit_ix_).init;
            
            x = obj.p2x(obj.s.xl,fit_init.p);
            x_lb = obj.p2x(obj.s.xl,fit_init.p_lb);
            x_ub = obj.p2x(obj.s.xl,fit_init.p_ub);
            
            if strcmpi(obj.opt.computeAlgo,'PS')
                
                minoptions = optimoptions(@patternsearch,'Display','iter',...
                    'MaxIter',obj.opt.MaxIter,'TolFun',obj.opt.TolX,'TolX',obj.opt.TolX,...
                    'MaxFunEvals',obj.opt.MaxIter*2,'StepTolerance',1e-3,...
                    'InitialMeshSize',2,'AccelerateMesh',obj.opt.ps_AccelerateMesh,...
                    'UseParallel',obj.opt.parallelsearch,'UseCompletePoll',obj.opt.parallelsearch);%,
                
                f = @(x) obj.opt.h_cost(x,fit_init.p);
                
                fprintf('Starting fit for %s\n',obj.subject);
                [x,Fps] = patternsearch(@(x)...
                    f(x)...
                    ,x,[],[],[],[],x_lb,x_ub,minoptions);
                
                if isfield(obj.opt,'dodoublefit')
                    dodoublefit = obj.opt.dodoublefit;
                else
                    dodoublefit = true;
                end
                
                if dodoublefit
                    fprintf('Repeat from optimum (should be quick...).\n');
                    [x,Fps] = patternsearch(@(x)...
                        f(x)...
                        ,x,[],[],[],[],x_lb,x_ub,minoptions);
                end
            else
                error('Invalid compute algo')
            end
            obj.fit(fit_ix_).options = minoptions;
            obj.fit(fit_ix_).p = obj.px2p(obj.s.xl,fit_init.p,x);
            
            [nll_app,aic_app,aicc_app,bic_app,ll_vec_app] = obj.opt.h_cost([],obj.fit(fit_ix_).p);
            obj.fit(fit_ix_).nll = nll_app;
            obj.fit(fit_ix_).aic = aic_app;
            obj.fit(fit_ix_).ll_vec = ll_vec_app;
            obj.fit(fit_ix_).n = sum(not(isnan(ll_vec_app)));
            obj.fit(fit_ix_).aicc = aicc_app;
            obj.fit(fit_ix_).bic = bic_app;
            
            %re-save this data with the fit
            obj.fit(fit_ix_).id_model = obj.id_model;
            obj.fit(fit_ix_).id_search = obj.id_search;
            obj.fit(fit_ix_).s = obj.s;
            obj.fit(fit_ix_).modelKey = obj.modelKey;
            obj.fit(fit_ix_).opt = obj.opt;
            obj.fit(fit_ix_).info = obj.info;
            %this is here, because we don't want to increment it if
            %something crashes before this function ends.
            obj.fit_ix = fit_ix_;
            
        end
        
        function ddm_mcmc(obj,varargin)
            d.mccount = 25e3;
            d.ThinChain = 5;
            d.doParallel = true;
            d.n_s = 25;
            d.BurnIn = 0.25;
            v = inputParser;
            addOptional(v,'mccount',d.mccount);
            addOptional(v,'ThinChain',d.ThinChain)
            addOptional(v,'doParallel',d.doParallel)
            addOptional(v,'n_s',d.n_s)
            addOptional(v,'BurnIn',d.BurnIn)
            parse(v,varargin{:})
            opt_ = v.Results;
            d = [];clear d;
            
            function ll = prior_likelihood(x,xl,h_priors)
                vec_xl = fieldnames(xl);
                ll = nan(1,length(vec_xl));
                for ix_vec_xl = 1:length(vec_xl)
                    x_index = xl.(vec_xl{ix_vec_xl});
                    x_value = x(x_index);
                    ll(x_index) = h_priors.(vec_xl{ix_vec_xl})(x_value);
                end
                ll = sum(log(ll));
            end
            
            function minit = ddm_mcmc_minit(obj,x,xl,n_s)
                %this is a bit ad-hoc, but I am getting the starting
                %chains proportional to how close they are to the optimal
                %found by patternsearch.
                % generate a bunch of candidate parameters
                fprintf('Running MCMC initialisation based on pattern search...\n');
                n_init_rand = 5e3;
                x_rand = nan(n_init_rand,obj.s.fit_n);
                for ix_draw_prior = 1:n_init_rand
                    x_rand(ix_draw_prior,:) = obj.p2x(xl,ddm_get_instance(obj, 'init_random'));
                end
                %get distance from ps optimum
                xd = x_rand-x;
                %normalise and convert to larger if closer to x
                xd = 1./rms(xd./std(xd,1),2);
                abc = cumsum(xd/sum(xd))-rand(1,n_s);
                %choose values proportional to their size
                abc(abc<0) = inf;
                [~,ix_x_rand] = min(abs(abc));
                minit = x_rand(ix_x_rand,:)';
                fprintf('Done.\n');
            end
            
            %use most uptodata p value
            p = obj.fit(obj.fit_ix).p;
            xl = obj.s.xl;
            h_priors = ddm_get_instance(obj, 'prior');
            x = obj.p2x(xl,p);
            
            minit = ddm_mcmc_minit(obj,x,xl,opt_.n_s);
            
            %negative because by default cost function returns the nll
            loglike = @(m) -obj.opt.h_cost(m,p);
            logprior = @(m) prior_likelihood(m, xl, h_priors);
            logPfuns = {@(m)logprior(m) @(m)loglike(m)};
            %             logPfuns{1}(minit(:,1))+logPfuns{2}(minit(:,1))
            fprintf('Starting mcmc...\n');
            [models,logp] = gwmcmc(minit,logPfuns,opt_.mccount,'Parallel',...
                opt_.doParallel,'BurnIn',opt_.BurnIn,'ThinChain',opt_.ThinChain);
            obj.mcmc(obj.fit_ix).models = models;
            obj.mcmc(obj.fit_ix).logp = logp;
            obj.mcmc(obj.fit_ix).opt = opt_;
            %should add a way to resume these I guess.
        end
        
        function get_data(obj)
            if not(isempty(obj.subject))
                data_ = readtable(obj.path_data);
                if isnumeric(data_.subject)
                    %if in csv subject id is numeric, then subject should
                    %be specified as 'sub01', 'sub02', etc.
                    data_.subject = arrayfun(@(x) sprintf('sub%02d',x),data_.subject,'UniformOutput',false);
                end
                case_subject = strcmpi(data_.subject,obj.subject);
                
                %                 if iscell(data_.choice)
                %                     error('Maybe you have a string in choice field?');
                %                 end
                data_ = data_(case_subject,:);
                if not(height(data_)>0),error('Bad subject spec?');end
                
                
                
                obj.data = data_;
            else
                error('No subject specified');
            end
        end
        
        
        function f_savepath = ddm_get_save_path(obj, f_path)
            if not(exist('f_path','var')==1),f_path = fullfile('sim',obj.modelclass);end
            if isempty(f_path),f_path = fullfile('sim',obj.modelclass);end
            f_name = sprintf('%s_%s_%s_%s%s.mat',...
                obj.subject,...
                ddm_def.debi_model(obj.id_model,'de','st'),...
                ddm_def.debi_model(obj.id_search,'de','st'),...
                obj.modelclass,...
                obj.extra_string);
            if not(exist(f_path,'dir')==7),mkdir(f_path);end
            f_savepath = fullfile(f_path,f_name);
        end
        function f_savepath = ddm_save(obj, f_path)
            if not(exist('f_path','var')==1),f_path = '';end
            f_savepath = ddm_get_save_path(obj, f_path);
            
            save(f_savepath,'obj');
        end
        
        function [f_savepath,obj_] = ddm_save_properties(obj, f_path)
            if not(exist('f_path','var')==1),f_path = '';end
            f_savepath = ddm_get_save_path(obj, f_path);
            f_savepath = strrep(f_savepath,'.mat','_properties.mat');
            
            prop_names = properties(obj);
            prop_names(strcmpi(prop_names,'objp')) = [];
            obj.objp = [];
            for ix_prop_names = 1:length(prop_names)
                obj.objp.(prop_names{ix_prop_names}) = obj.(prop_names{ix_prop_names});
            end
            obj_ = obj.objp;
            save(f_savepath,'obj_');
        end
        
        function [t_ups,t_dow,p_ups,p_dow] = ddm_data_draw(obj,p,N)
            if contains(func2str(obj.ddm_pdf),'ddm_prt_ana')
                error(['Need a full pdf generator from which to draw.', ...
                    'Change the pdf function handle.',...
                    'e.g.: obj.ddm_pdf = @(a,b) obj.ddm_pdf_ana(a,b);']);
            end
            lt  = 0:obj.s.dt:obj.s.T-obj.s.dt;
            
            [~,~,t,cdf_dow,cdf_ups] = obj.ddm_pdf(p,lt);
            [t_ups,t_dow,p_ups,p_dow] = obj.ddm_pdf2rt(cdf_ups,cdf_dow,t,N);
        end
        
        
        
        function outputArg = ddm_get_instance(obj, deftype)
            [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ...
                obj.ddm_def_instance;
            
            for ix_modelkey_var = 1:length(modelkey_var)
                modelkey_rev.(modelkey_var{ix_modelkey_var}) = ix_modelkey_var;
            end
            
            if strcmpi(deftype,'keyf')
                outputArg = modelkey_var;
            elseif strcmpi(deftype,'keyr')
                outputArg = modelkey_rev;
            elseif strcmpi(deftype,'init_random')
                outputArg = pran_;
            elseif strcmpi(deftype,'init_default')
                outputArg = pdef_;
            elseif strcmpi(deftype,'lbound')
                outputArg = plbound_;
            elseif strcmpi(deftype,'ubound')
                outputArg = pubound_;
            elseif strcmpi(deftype,'prior')
                outputArg = prior_;
            else
                error('bad deftype')
            end
        end
        
        
        
        function p_mat = ddm_cost_add_stim_dependencies(obj,p_mat)
            %             p_mat.c = obj.data.stim_conflict;
        end
        
        function h_f = aux_plot(obj,varargin)
            % Plot function for pdf/data - a lot like cost function
            % At the moment this function would deal terribly
            % with individual trial modulation data - (each trial would
            % spawn a subplot).
            % Example use:
            % %             for ix_sub = 1:length(sub_cell)
            % %                 c = cmap(issz(ix_sub)+1,:);
            % %                 es = sprintf('_isssz%d',issz(ix_sub));
            % %                 skipcases.difficulty = [-2,-1,+1,+2];
            % %                 h_f = sr_cell{ix_sub}.aux_plot('color',c,'es',es,'plotsinglerow',true,'skipcases',skipcases);
            % %                 
            % %                 drawnow;
            % %                 printForPub(gcf,h_f.Name,...
            % %                     'fformat',fig.fformat,'physicalSizeCM',[25,5],'savedir',fig.savedir,'doPrint',fig.doPrint);
            % %             end
            %
            addpath('ext');%want to use numSubplots function - not sure if this points to the right folder though
            %
            d.color = [0 0 0];
            d.es = '';%extra string in title name
            d.plotsinglerow = false;
            d.p = obj.fit(end).p;
            d.skipcases = struct;
            %
            v = inputParser;
            addOptional(v,'color',d.color)
            addOptional(v,'es',d.es)
            addOptional(v,'plotsinglerow',d.plotsinglerow)
            addOptional(v,'p',d.p)
            addOptional(v,'skipcases',d.skipcases)
            parse(v,varargin{:})
            v = v.Results;
            d = [];clear d;
            %
            p = v.p;v = rmfield(v,'p');
            %
            p_mat = struct2table(repmat(p,height(obj.data),1));
            
            %Add trial-trial stimulus dependencies (e.g. coherence, conflict etc.)
            p_mat = obj.ddm_cost_add_stim_dependencies(p_mat);
            
            p_mat_unique = unique(p_mat);
            p_mat_array = table2array(p_mat);
            
            case_nan = isnan(obj.data.rt)|isnan(obj.data.choice);
            case_right = obj.data.choice;
            %this is ugly - setting nan to zero so that I can use them as
            %logicals.case_nnan is used to exclude them later though so ok.
            case_right(case_nan) = 0;
            case_right = logical(case_right);
            case_wrong = not(case_right);
            
            %skipcases!
            fn = fieldnames(v.skipcases);
            for ix_fn = 1:length(fn)
                p_mat_unique(ismember(p_mat_unique.(fn{ix_fn}),v.skipcases.(fn{ix_fn})),:) = [];
            end
            
            if v.plotsinglerow
                ns = [1,height(p_mat_unique)];
            else
                [ns,~] = numSubplots(height(p_mat_unique));
            end
            lt  = 0:obj.s.dt:obj.s.T-obj.s.dt;
            
            h_f = figure;clf;
            h_f.Name = sprintf('%s_%s%s',obj.modelclass,obj.subject,v.es);
            %find which variable is changing
            ix_change = find(arrayfun(@(ix) height(unique(p_mat_unique(:,ix))),1:size(p_mat_unique,2))>1);
            if length(ix_change)>1,warning('Not equipped to deal with multiple changes');
            else
            name_change = p_mat_unique.Properties.VariableNames{ix_change};
            name_change_label = name_change;
            if length(name_change_label)>4,name_change_label = [name_change_label(1:4) '.'];end
            end
            
            
            for ix_p_config = 1:height(p_mat_unique)
                px = table2struct(p_mat_unique(ix_p_config,:));
                px_array = table2array(p_mat_unique(ix_p_config,:));
                case_config = all(p_mat_array==px_array,2);
                t_cr = obj.data.rt(case_right&case_config&not(case_nan));
                t_cw = obj.data.rt(case_wrong&case_config&not(case_nan));
                
                title_ = sprintf('%s:\n%0.1f',name_change_label,px.(name_change));
                %Get simulated pdf
                if contains(func2str(obj.ddm_pdf),'ddm_prt_ana')
                    [pdf_dow,pdf_ups,rt,~,~] = obj.ddm_pdf_ana(px,lt);
                else
                    [pdf_dow,pdf_ups,rt,~,~] = obj.ddm_pdf(px,lt);
                end
                
                %Get data pdf
                bw = 0.065;
                fr = length(t_cr)/(length(t_cr)+length(t_cw));
                if length(t_cr)>1
                    pdf_ups_dat = ksdensity(t_cr,rt,'Bandwidth',bw)*fr;
                else
                    pdf_ups_dat = zeros(size(rt));
                end
                if length(t_cw)>1
                    pdf_dow_dat = ksdensity(t_cw,rt,'Bandwidth',bw)*(1-fr);
                else
                    pdf_dow_dat = zeros(size(rt));
                end
                
                %ugly (remove zero before flip)
                rt0 = rt==0;rt(rt0) = [];
                pdf_dow(rt0) = [];pdf_ups(rt0) = [];
                pdf_dow_dat(rt0) = [];pdf_ups_dat(rt0) = [];
                
                rte = [-fliplr(rt),rt];
                pdf = [fliplr(pdf_dow_dat),pdf_ups_dat];
                pdf_dat = [fliplr(pdf_dow),pdf_ups];
                
                if not(length(unique(rte))==2*length(rt))
                    error('Zero in the middle? Need to fix.');
                end
                
                subplot(ns(1),ns(2),ix_p_config);cla;hold on;
                h1 = plot(rte,pdf_dat,'Color',v.color,'lineWidth',2);
                h2 = plot(rte,pdf,'--','Color',0.6*ones(1,3),'lineWidth',2);
                title(title_);
                axis tight;
            end
            legend([h1,h2],{'Data','Sim'},'Location','NorthEast');
            xlabel('RT (s)')
        end
        
        function [nll_app,aic_app,aicc_app,bic_app,ll_vec] = ddm_cost_pdf_nll(obj,x,p)
            if not(isempty(x))
                p = obj.px2p(obj.s.xl,p,x);
            end
            p_RT_and_accuracy = nan(height(obj.data),1);
            
            p_mat = struct2table(repmat(p,height(obj.data),1));
            
            %Add trial-trial stimulus dependencies (e.g. coherence, conflict etc.)
            p_mat = obj.ddm_cost_add_stim_dependencies(p_mat);
            
            p_mat_unique = unique(p_mat);
            p_mat_array = table2array(p_mat);
            
            case_nan = isnan(obj.data.rt)|isnan(obj.data.choice);
            case_right = obj.data.choice;
            %this is ugly - setting nan to zero so that I can use them as
            %logicals.case_nnan is used to exclude them later though so ok.
            case_right(case_nan) = 0;
            case_right = logical(case_right);
            case_wrong = not(case_right);
            lt  = 0:obj.s.dt:obj.s.T-obj.s.dt;
            for ix_p_config = 1:height(p_mat_unique)
                px = table2struct(p_mat_unique(ix_p_config,:));
                px_array = table2array(p_mat_unique(ix_p_config,:));
                case_config = all(p_mat_array==px_array,2);
                t_cr = obj.data.rt(case_right&case_config&not(case_nan));
                t_cw = obj.data.rt(case_wrong&case_config&not(case_nan));
                
                if contains(func2str(obj.ddm_pdf),'ddm_prt_ana')
                    
                    [pRT_g_cr, p_cr] = obj.ddm_pdf(px,+t_cr);
                    [pRT_g_cw, ~] = obj.ddm_pdf(px,-t_cw);
                    pRT_g_cr_x_p_cr = pRT_g_cr*p_cr;
                    pRT_g_cw_x_p_cw = pRT_g_cw*(1-p_cr);
                    
                    
                else
                    ix_cr = round(t_cr/obj.s.dt);
                    ix_cw = round(t_cw/obj.s.dt);
                    [pdf_cw,pdf_cr,~,cdf_cw,cdf_cr] = obj.ddm_pdf(px,lt);
                    p_cr = cdf_cr(end);
                    p_cw = cdf_cw(end);
                    pRT_g_cr_x_p_cr = pdf_cr(ix_cr)'*p_cr;
                    pRT_g_cw_x_p_cw = pdf_cw(ix_cw)'*p_cw;
                    
                end
                
                p_RT_and_accuracy(case_right&case_config&not(case_nan)) = pRT_g_cr_x_p_cr;
                p_RT_and_accuracy(case_wrong&case_config&not(case_nan)) = pRT_g_cw_x_p_cw;
            end
            p_RT_and_accuracy(p_RT_and_accuracy == 0) = 1e-32;%not great
            
            
            ll_vec = log(p_RT_and_accuracy);
            ll_app = nansum(ll_vec);%nan are values that are not in condition
            %     if isnan(ll_app)||isinf(ll_app),error('non scallr ll');end
            if isnan(ll_app),ll_app=-inf;end
            k = obj.s.fit_n+1;%number of free params + 1 for fixed noise
            n = sum(not(isnan(p_RT_and_accuracy)));
            %prob an issue here with thr fact that some trials are kicked out..
            bic_app = log(n)*k-2*ll_app;
            aic_app = 2*k-2*ll_app;
            aicc_app = aic_app + (2*(k^2) + 2*k)/(n-k-1);
            nll_app = -ll_app;
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance(obj)
            ix = 1;
            
            p_ = 's';
            modelkey_var{ix} = p_;ix = ix+1;
            pran_.(p_) = 1;
            pdef_.(p_) = 1;
            plbound_.(p_) = 1;
            pubound_.(p_) = 1;
            prior_.(p_) = @(x) 1;
            
            p_ = 'z';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_alpha = 5;
            g_beta = 5;
            pd_hn = makedist('beta','a',g_alpha,'b',g_beta);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.5;
            plbound_.(p_) = 0;
            pubound_.(p_) = 1;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'v';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 3;g_sd = 1;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_sd);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 0.0;
            plbound_.(p_) = -20;
            pubound_.(p_) = +20;
            prior_.(p_) = @(x) gampdf(x,A_shape,B_scale);
            
            p_ = 'a';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 1;g_sd = 0.25;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_sd);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 1.0;
            plbound_.(p_) = 0.1;
            pubound_.(p_) = 7.5;
            prior_.(p_) = @(x) gampdf(x,A_shape,B_scale);
            
            p_ = 't';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 0.3;g_sd = 0.075;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_sd);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 0.25;
            plbound_.(p_) = 0.1;
            pubound_.(p_) = 1.5;
            prior_.(p_) = @(x) gampdf(x,A_shape,B_scale);
            
            p_ = 'st';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_lo = 0;
            g_up = 0.25;
            pran_.(p_) = unifrnd(g_lo,g_up);
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 1.0;
            prior_.(p_) = @(x) unifpdf(x,g_lo,g_up);
            
            p_ = 'sv';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 2;
            pd_hn = makedist('HalfNormal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 10;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'sz';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_alpha = 1;
            g_beta = 8;
            pd_hn = makedist('beta','a',g_alpha,'b',g_beta);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 0.75;%the real bound is 1... but that's a bit much
            %esp since generally looks like it's more like 0.01 if anything
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
        end
    end
    methods (Static)
        
        function  [pdf_,p_cr] = ddm_prt_ana(p,rt)
            err = 1e-8;
            %it's a surprise to me that other terms don't matter here...
            p_cr = ddm_def.hddm_prob_ub(p.v,p.a,p.z);
            h_pdf = @(x) ddm_def.hddm_pdf_full(x,p.v,p.sv,p.a,p.z,p.sz,p.t,p.st,err);
            pdf_ = arrayfun(@(x) h_pdf(x),+rt);
            
        end
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_ana(p,rt)
            if any(rt<0),error(['-rt is only usable for ddm_prt_ana.\n' ...
                    'In normal pdf functions you ask for +rt, but you get '...
                    'out pdf for both choices. Maybe this should be changed '...
                    'in future. Also this check probably takes time since '...
                    'rt is often large.']);
            end
            
            err = 1e-8;
            
            h_pdf = @(x) ddm_def.hddm_pdf_full(x,p.v,p.sv,p.a,p.z,p.sz,p.t,p.st,err);
            pdf_ups = arrayfun(@(x) h_pdf(x),+rt);
            pdf_dow = arrayfun(@(x) h_pdf(x),-rt);
            
            
            cdf_ups = cumtrapz(rt, pdf_ups);
            cdf_dow = cumtrapz(rt, pdf_dow);
            p_ups = ddm_def.hddm_prob_ub(p.v,p.a,p.z);
            %           Rescale the cdf to avoid numerical issues with the integration
            %           assumes that max(rt) really has made the tail to go zero.
            cdf_ups = (cdf_ups/cdf_ups(end))*p_ups;
            cdf_dow = (cdf_dow/cdf_dow(end))*(1-p_ups);
            
            
        end
        
        function [cdf_ups,cdf_dow] = cdf_t_st(cdf_ups,cdf_dow,rt,t,st,dt)
            td0_vec = find((rt>(t-st/2))&(rt<(t+st/2)));
            if length(td0_vec)<=1
                ix_t0_shift = round((t)/dt);
                cdf_ups = circshift(cdf_ups,ix_t0_shift);
                cdf_ups(1:ix_t0_shift) = 0;
                cdf_dow = circshift(cdf_dow,ix_t0_shift);
                cdf_dow(1:ix_t0_shift) = 0;
            else
                N_vec_c_td = length(td0_vec);
                cdf_ups_mats = repmat(cdf_ups,N_vec_c_td,1)/N_vec_c_td;
                cdf_dow_mats = repmat(cdf_dow,N_vec_c_td,1)/N_vec_c_td;
                
                for ix_td0_vec = 1:length(td0_vec)
                    ix_t0_shift = td0_vec(ix_td0_vec);
                    cdf_ups_mats(ix_td0_vec,:) = circshift(cdf_ups_mats(ix_td0_vec,:),ix_t0_shift);
                    cdf_ups_mats(ix_td0_vec,1:ix_t0_shift) = 0;
                    cdf_dow_mats(ix_td0_vec,:) = circshift(cdf_dow_mats(ix_td0_vec,:),ix_t0_shift);
                    cdf_dow_mats(ix_td0_vec,1:ix_t0_shift) = 0;
                end
                cdf_ups = sum(cdf_ups_mats,1);
                cdf_dow = sum(cdf_dow_mats,1);
                
            end
        end
        function [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_bru(p,lt,N_its)
            
            dt = lt(2)-lt(1);
            T = lt(end)+dt;
            maxIterations = floor(T/dt);
            %             if (length(lt)-1)~=maxIterations,error('Messed up time code');end
            if (length(lt))~=maxIterations,error('Messed up time code');end
            %%
            x0 = p.z*p.a;
            x0_trial = (rand(N_its,1)-0.5)*(p.sz) + x0;%n.b. this one is uniform.
            x_noise = randn(N_its,maxIterations)*(p.s)*sqrt(dt);
            x_drift = repmat((p.v + p.sv*randn(N_its,1))*dt,1,maxIterations);
            %%
            x = nan(N_its,maxIterations);
            x(:,1) = x0_trial;
            for ix_t = 2:maxIterations
                %need to check this before it gets used:
                x(:,ix_t) = x(:,ix_t-1) + x_noise(:,ix_t) + x_drift(:,ix_t);
            end
            %%
            x_threshold = +p.a;
            [~,x_upper]=sort(x>x_threshold,2,'descend');% should the sign here be the same????
            x_upper=x_upper(:,1);
            x_upper(x_upper==1) = nan;
            
            x_threshold = 0;
            [~,x_lower]=sort(x<x_threshold,2,'descend');% should the sign here be the same????
            x_lower=x_lower(:,1);
            x_lower(x_lower==1) = nan;
            
            x_bounds = [x_upper,x_lower];
            x_bounds = gather(x_bounds);
            [ix_time,vec_correct] = nanmin(x_bounds,[],2);
            ix_time(isnan(ix_time)) = length(lt);%not sure if this is best way to deal with problem...
            vec_rt = lt(ix_time)';
            % don't think the time shift can change winner/loser so do it here.
            
            td_tot = p.t + (rand(N_its,1)-0.5)*p.st;
            vec_rt = vec_rt + td_tot;
            
            vec_correct = vec_correct==1;
            
            cdf_ups = ksdensity(vec_rt(vec_correct),lt,'Support','Positive','function','cdf');
            cdf_dow = ksdensity(vec_rt(not(vec_correct)),lt,'Support','Positive','function','cdf');
            p_ups = sum(vec_correct)/length(vec_correct);
            p_dow = 1-p_ups;
            %this is already handled at the moment.
            
            cdf_ups = cdf_ups*p_ups;
            cdf_dow = cdf_dow*p_dow;
            %
            pdf_ups = diff(cdf_ups)/dt;
            pdf_dow = diff(cdf_dow)/dt;
            rt = (lt(1:end-1)+lt(2:end))*0.5;
        end
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_trm(p,lt,dx)
            %
            dt = lt(2)-lt(1);
            f = 5;%increasing this helps (as does decreasing dx)
            xmax = p.a + f*sqrt(dt)*p.s;
            xmin = 0 - f*sqrt(dt)*p.s;
            
            xz = [xmax:-dx:xmin]';
            xvm_probe = repmat(xz',length(xz),1)';
            xvm_prev = repmat(xz',length(xz),1);
            
            za = (p.z*p.a);
            [~,zeroStateIx] = min(abs(xz-za));
            if (p.sz <= 1e-3)
                x0 = zeros(length(xz),1);
                x0(zeroStateIx,1) = 1;
            else
                x0 = unifpdf(xz,xz(zeroStateIx)-p.sz/2,xz(zeroStateIx)+p.sz/2);
                x0 = x0/sum(x0);
            end
            %
            sig = p.s*sqrt(dt);
            N_t = length(lt);
            if p.sv<1e-3%arb. threshold used in hddm core
                p.sv = 0;
                N_sv = 1;
                vec_sv = 0;
                p_sv = 1;
            else
                N_sv = 11;
                vec_sv = p.sv*linspace(-2.5,2.5,N_sv);
                p_sv = normpdf(vec_sv,0,p.sv);
                p_sv = p_sv/sum(p_sv);
            end
            
            % CORE
            cdf_ups = nan(N_sv,N_t);
            cdf_dow = nan(N_sv,N_t);
            for ix_sv = 1:N_sv
                v = (p.v+vec_sv(ix_sv));
                pMat = zeros(length(xz),N_t);
                
                zn = x0;%used to have a multiplication with v here (not sure why)
                xz_ups = xz>p.a;
                xz_dow = xz<0;
                e_ups = eye(sum(xz_ups));
                e_dow = eye(sum(xz_dow));
                len_e_ups = length(e_ups);
                len_e_dow = length(e_dow);
                
                pMat(:,1) = zn;
                for ix_t = 2:N_t
                    xvm_expect = xvm_prev + v*dt;
                    A = (1/sqrt(2*pi*(sig^2)))*exp(-((xvm_probe - xvm_expect).^2)/(2*(sig^2)));
                    An = A./repmat(sum(A,1),length(xz),1);
                    An(isnan(An))=0;
                    
                    An(:,xz>p.a) = 0;
                    An(1:len_e_ups,1:len_e_ups) = e_ups;
                    An(:,xz<0) = 0;
                    An(end-len_e_dow+1:end,end-len_e_dow+1:end) = e_dow;
                    
                    zn = An*zn;
                    pMat(:,ix_t) = zn;
                end
                cdf_ups(ix_sv,:) = sum(pMat(xz_ups,:));
                cdf_dow(ix_sv,:) = sum(pMat(xz_dow,:));
            end
            
            cdf_dow = p_sv*cdf_dow;
            cdf_ups = p_sv*cdf_ups;
            %
            [cdf_ups,cdf_dow] = ddm_def.cdf_t_st(cdf_ups,cdf_dow,lt,p.t,p.st,dt);
            
            cdf_ups_end = cdf_ups(end);
            cdf_dow_end = cdf_dow(end);
            cdf_ups = cdf_ups/(cdf_ups_end+cdf_dow_end);
            cdf_dow = cdf_dow/(cdf_ups_end+cdf_dow_end);
            %
            pdf_ups = diff(cdf_ups)/dt;
            pdf_dow = diff(cdf_dow)/dt;
            rt = 0.5*(lt(1:end-1)+lt(2:end));
        end
        
        function [t_ups,t_dow,p_ups,p_dow] = ddm_pdf2rt(cdf_ups,cdf_dow,t_math,N,varargin)
            d.balanced = true;
            %%
            v = inputParser;
            addOptional(v,'balanced',d.balanced)
            parse(v,varargin{:})
            v = v.Results;
            d = [];clear d;
            %%
            p_ups = (cdf_ups(end)/(cdf_ups(end)+cdf_dow(end)));
            p_dow = 1 - p_ups;
            if v.balanced
                N_ups = round(N*p_ups);
                N_dow = N - N_ups;
            else
                N_ups = N;
                N_dow = N;
            end
            
            r_ups = rand(N_ups,1)*p_ups;
            r_dow = rand(N_dow,1)*(1-p_ups);
            [~,ix_ups] = min(abs(cdf_ups-r_ups),[],2);
            [~,ix_dow] = min(abs(cdf_dow-r_dow),[],2);
            t_ups = t_math(ix_ups);
            t_dow = t_math(ix_dow);
        end
        
        function [A_shape,B_scale] = gamma_convert(g_mea,g_std)
            alpha_shape = (g_mea^2)/(g_std^2);
            beta_rate = (g_mea)/(g_std^2);
            A_shape = alpha_shape;
            B_scale = 1/beta_rate;
        end
        
        
        function x = p2x(xl,p)
            vec_xl = fieldnames(xl);
            x = nan(1,length(vec_xl));
            for ix_vec_xl = 1:length(vec_xl)
                x_index = xl.(vec_xl{ix_vec_xl});
                x_value = p.(vec_xl{ix_vec_xl});
                x(x_index) = x_value;
            end
        end
        
        function p = px2p(xl,p,x)
            vec_xl = fieldnames(xl);
            for ix_vec_xl = 1:length(vec_xl)
                x_index = xl.(vec_xl{ix_vec_xl});
                x_value = x(x_index);
                p.(vec_xl{ix_vec_xl}) = x_value;
            end
        end
        
        function op = debi_model(ip, ip_type, op_req)
            
            max_n_param = 39;
            
            if strcmpi(ip_type,'de')
                if ip>((2^(max_n_param+1))-1),error('bad de debi use');end
            elseif strcmpi(ip_type,'bi')
                if length(ip)>max_n_param,error('bad bi debi use');end
            end
            %                         if strcmpi(op_req,'de')
            %             elseif strcmpi(op_req,'bi')
            %             end
            
            %core_char_length = max decimal length of binary number
            %             core_char_length = length(num2str(bi2de(ones(1,max_n_param),'left-msb')));
            if strcmpi(ip_type,'de')&&strcmpi(op_req,'bi')
                op = de2bi(ip,max_n_param,'left-msb');
            elseif strcmpi(ip_type,'bi')&&strcmpi(op_req,'de')
                op = bi2de(ip,'left-msb');
            elseif strcmpi(ip_type,'bi')&&strcmpi(op_req,'st')
                tempX = bi2de(ip,'left-msb');
                op = ['x' sprintf('%010d',tempX) 'x'];
            elseif strcmpi(ip_type,'de')&&strcmpi(op_req,'st')
                op = ['x' sprintf('%010d',ip) 'x'];
            end
            
            if strcmpi(op_req,'de')
                if op>((2^(max_n_param+1))-1),error('bad de debi use');end
            elseif strcmpi(op_req,'bi')
                if length(op)>max_n_param,error('bad bi debi use');end
            end
        end
        
        %
        function p = ftt_01w(tt,w,err)
            %Translated from HDDM core/N.F paper
            if (pi*tt*err)<1 % if error threshold is set low enough
                kl=sqrt(-2*log(pi*tt*err)/(pi^2*tt)); % bound
                kl=max(kl,1/(pi*sqrt(tt))); % ensure boundary conditions met
            else % if error threshold set too high
                kl=1/(pi*sqrt(tt)); % set to boundary condition
            end
            % calculate number of terms needed for small t
            if (2*sqrt(2*pi*tt)*err)<1 % if error threshold is set low enough
                ks=2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)); % bound
                ks=max(ks,sqrt(tt)+1); % ensure boundary conditions are met
            else % if error threshold was set too high
                ks=2; % minimal kappa for that case
            end
            
            % compute f(tt|0,1,w)
            p=0; %initialize density
            if ks<kl % if small t is better (i.e., lambda<0)
                K=ceil(ks); % round to smallest integer meeting error
                vec_k = -floor((K-1)/2):ceil((K-1)/2);
                for k = vec_k
                    p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
                end
                p=p/sqrt(2*pi*(tt^3)); % add con_stant term
                
            else % if large t is better...
                K= ceil(kl); % round to smallest integer meeting error
                
                for k=1:K
                    p=p+k*exp(-((k^2))*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
                end
                p=p*pi; % add con_stant term
            end
        end
        
        function op = hddm_prob_ub(v,a,z)
            %"""Probability of hitting upper boundary."""
            op = (exp(-2*a*z*v) - 1) / (exp(-2*a*v) - 1);
        end
        
        
        function op = hddm_pdf(x,v,a,w,err)
            %Translated from HDDM core
            %"""Compute the likelihood of the drift diffusion model f(t|v,a,z) using the method
            %and implementation of Navarro & Fuss, 2009.
            %"""
            if x <= 0
                op = 0;
            else
                
                tt = x/(a^2);% # use normalized time
                p = ddm_def.ftt_01w(tt, w, err);% #get f(t|0,1,w)
                
                %% convert to f(t|v,a,w)
                op =  p*exp(-v*a*w -(v^2)*x/2)/(a^2);
            end
        end
        
        function op = hddm_pdf_sv(x,v,sv,a,z,err)
            %Translated from HDDM core
            %"""Compute the likelihood of the drift diffusion model f(t|v,a,z,sv) using the method
            %and implementation of Navarro & Fuss, 2009.
            %sv is the std of the drift rate
            %"""
            if x <= 0
                op = 0;
            elseif sv==0
                op =  ddm_def.hddm_pdf(x, v, a, z, err);
            else
                tt = x/(a^2); %% use normalized time
                p  = ddm_def.ftt_01w(tt, z, err); %%get f(t|0,1,w)
                
                %% convert to f(t|v,a,w)
                op = exp(log(p) + ((a*z*sv)^2 - 2*a*v*z - (v^2)*x)/(2*(sv^2)*x+2))/sqrt((sv^2)*x+1)/(a^2);
            end
        end
        
        
        
        function op = hddm_pdf_full(x,v,sv,a,z,sz,t,st,err)
            
            %Translated from HDDM core
            n_st = 15;
            n_sz = 15;
            use_adaptive = 0;
            %             simps_err = 1e-8;
            %"""full pdf"""
            %% Check if parpameters are valid (also don't compute if less than t)
            if (z<0) || (z>1) || (a<0) || (t<0) || (st<0) || (sv<0) || (sz<0) || (sz>1) ||...
                    ((abs(x)-(t-st/2.))<0) || (z+sz/2.>1) || (z-sz/2.<0) || (t-st/2.<0)
                op = 0;
                return;
            end
            
            %% transform x,v,z if x is upper bound response
            if x > 0
                v = -v;
                z = 1.-z;
            end
            x = abs(x);
            
            if st<1e-3
                st = 0;
            end
            if sz <1e-3
                sz = 0;
            end
            
            if (sz==0)
                if (st==0) %%sv=$,sz=0,st=0
                    op = ddm_def.hddm_pdf_sv(x - t, v, sv, a, z, err);
                else      %%sv=$,sz=0,st=$
                    if use_adaptive>0
                        error('Not copied over');
                    else
                        op = ddm_def.hddm_simpson_1D(x, v, sv, a, z, t, err, z, z, 0, t-st/2., t+st/2., n_st);
                    end
                end
            else %%sz=$
                if (st==0) %%sv=0,sz=$,st=0
                    if use_adaptive
                        error('Not copied over');
                    else
                        op = ddm_def.hddm_simpson_1D(x, v, sv, a, z, t, err, z-sz/2., z+sz/2., n_sz, t, t , 0);
                    end
                else      %%sv=0,sz=$,st=$
                    if use_adaptive
                        error('Not copied over');
                    else
                        op = ddm_def.hddm_simpson_2D(x, v, sv, a, z, t, err, z-sz/2., z+sz/2., n_sz, t-st/2., t+st/2., n_st);
                    end
                end
            end
        end
        
        
        
        function op = hddm_simpson_1D(x, v, sv, a, z, t, err, lb_z, ub_z, n_sz, lb_t, ub_t, n_st)
            %Translated from HDDM core
            %assert ((n_sz&1)==0 and (n_st&1)==0), "n_st and n_sz have to be even"
            %assert ((ub_t-lb_t)*(ub_z-lb_z)==0 and (n_sz*n_st)==0), "the function is defined for 1D-integration only"
            
            n = max(n_st, n_sz);
            if n_st==0 %integration over z
                hz = (ub_z-lb_z)/n;
                ht = 0;
                lb_t = t;
                ub_t = t;
            else %integration over t
                hz = 0;
                ht = (ub_t-lb_t)/n;
                lb_z = z;
                ub_z = z;
            end
            S = ddm_def.hddm_pdf_sv(x - lb_t, v, sv, a, lb_z, err);
            
            
            for i = 1:n
                z_tag = lb_z + hz * i;
                t_tag = lb_t + ht * i;
                y = ddm_def.hddm_pdf_sv(x - t_tag, v, sv, a, z_tag, err);
                if mod(i,2)==1 %check if i is odd
                    S = S  + (4 * y);
                else
                    S = S + (2 * y);
                end
            end
            
            S = S - y; %the last term should be f(b) and not 2*f(b) so we subtract y
            
            S = S / ((ub_t-lb_t)+(ub_z-lb_z)); %the right function if pdf_sv()/sz or pdf_sv()/st
            
            op = ((ht+hz) * S / 3);
        end
        
        function op = hddm_simpson_2D(x, v, sv, a, z, t, err, lb_z, ub_z, n_sz, lb_t, ub_t, n_st)
            %Translated from HDDM core
            %assert ((n_sz&1)==0 and (n_st&1)==0), "n_st and n_sz have to be even"
            %assert ((ub_t-lb_t)*(ub_z-lb_z)>0 and (n_sz*n_st)>0), "the function is defined for 2D-integration only, lb_t: %f, ub_t %f, lb_z %f, ub_z %f, n_sz: %d, n_st %d" % (lb_t, ub_t, lb_z, ub_z, n_sz, n_st)
            
            ht = (ub_t-lb_t)/n_st;
            S = ddm_def.hddm_simpson_1D(x, v, sv, a, z, lb_t, err, lb_z, ub_z, n_sz, 0, 0, 0);
            
            for i_t  = 1:n_st
                t_tag = lb_t + ht * i_t;
                y = ddm_def.hddm_simpson_1D(x, v, sv, a, z, t_tag, err, lb_z, ub_z, n_sz, 0, 0, 0);
                if mod(i_t,2)==1 %check if i is odd
                    S = S + (4 * y);
                else
                    S = S + (2 * y);
                end
            end
            S = S - y; %the last term should be f(b) and not 2*f(b) so we subtract y
            S = S/ (ub_t-lb_t);
            
            op = (ht * S / 3);
        end
    end
end





