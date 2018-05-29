classdef ddm_def < matlab.mixin.Copyable%instead of handle
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        subject = '';
        modelclass = 'base';
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
    end
    
    methods
        function obj = ddm_def
            %light initialisation so functions can be used easily
            obj.modelclass = '';
            obj.modelKey = ddm_get_instance(obj, 'keyf');
            obj.info.version = sprintf('0.0.1');
            obj.info.date = datetime;
        end
        
        function ddm_init(obj, id_model,id_search)
            %Set the simulation settings (s), and the optimisation settings
            %(opt). id_model and if_fit (which are decimals,
            % which correspond to binary vectors referencing the
            % parameters) get set.
            % additionally id_search gets translated into xl which holds the
            % fit parameters as strings in a cell.
            obj.id_model = id_model;
            obj.id_search = id_search;
            
            if not(isnumeric(id_model)&(numel(id_model)==1)),error('Bad id_model');end
            if not(isnumeric(id_search)&(numel(id_search)==1)),error('Bad id_model');end
            
            obj.s.fit_n = sum(obj.debi_model(obj.id_search,'de','bi'));
            obj.s.minAlgo = 'nll';
            obj.s.reinit = false;
            obj.s.dt = 1e-3;
            obj.s.ddx = 100;
            obj.s.T = 5;
            obj.s.inittype = 'random';
            obj.s.path_data = '';
            
            id_search_index  = find(obj.debi_model(obj.id_search,'de','bi'));
            
            % Set the parameters over which to optimise from modelType spec
            co = 1;
            for ix_id_search_index = 1:length(id_search_index)
                parameter_string = obj.modelKey{id_search_index(ix_id_search_index)};
                obj.s.xl.(parameter_string) = co;co = co+1;
            end
            
            obj.opt.TolX = 1e-5;
            obj.opt.MaxIter = 1000;
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
            
            %             if not(length(obj.fit)==obj.fit_ix)
            %                 error('We already have an initialised fit here');
            %             end
            %             obj.fit_ix = obj.fit_ix + 1;
            
            if fit_ix_==1
                %p is defined from a full set of defaults/random init
                init_p_full = obj.ddm_get_instance(['init_' obj.s.inittype]);
                % and now we need to reduce it to match the specific model
                % configuration we are testing.
                id_model_index = find(obj.debi_model(obj.id_model,'de','bi'));
                for ix_parameter_cell = 1:length(obj.modelKey)
                    parameter_string = obj.modelKey{ix_parameter_cell};
                    if any(strcmp(parameter_string,obj.modelKey(id_model_index)))
                        init_p_reduced.(parameter_string) = init_p_full.(parameter_string);
                    else
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
            fit_init.nll = obj.opt.h_cost([],fit_init.p);
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
                [x,Fps] = patternsearch(@(x)...
                    f(x)...
                    ,x,[],[],[],[],x_lb,x_ub,minoptions);
            else
                error('Invalid compute algo')
            end
            obj.fit(fit_ix_).options = minoptions;
            obj.fit(fit_ix_).p = obj.px2p(obj.s.xl,fit_init.p,x);
            
            [nll_app,aic_app,aicc_app,bic_app] = obj.opt.h_cost([],obj.fit(fit_ix_).p);
            obj.fit(fit_ix_).nll = nll_app;
            obj.fit(fit_ix_).aic = aic_app;
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
            [models,logp]=gwmcmc(minit,logPfuns,opt_.mccount,'Parallel',opt_.doParallel,'BurnIn',opt_.BurnIn,'ThinChain',opt_.ThinChain);
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
                
                if iscell(data_.choice)
                    error('Maybe you have a string in choice field?');
                end
                data_ = data_(case_subject,:);
                if not(height(data_)>0),error('Bad subject spec?');end
                
                
                
                obj.data = data_;
            end
        end
        
        function f_savepath = ddm_save(obj,f_path)
            if not(exist('f_path','var')==1),f_path = fullfile('sim',obj.modelclass);end
            f_name = sprintf('%s_%s_%s_%s%s.mat',...
                obj.subject,...
                obj.debi_model(obj.id_model,'de','st'),...
                obj.debi_model(obj.id_search,'de','st'),...
                obj.modelclass,...
                obj.extra_string);
            if not(exist(f_path,'dir')==7),mkdir(f_path);end
            f_savepath = fullfile(f_path,f_name);
            save(f_savepath,'obj');
        end
        
        function [t_ups,t_dow,p_ups,p_dow] = ddm_data_draw(obj,p,N)
            [~,~,t,cdf_dow,cdf_ups] = obj.ddm_pdf(p,obj.s.dt,obj.s.T,obj.s.ddx);
            [t_ups,t_dow,p_ups,p_dow] = obj.ddm_pdf2rt(cdf_ups,cdf_dow,t,N);
        end
        
        function op = debi_model(obj, ip, ip_type, op_req)
            if strcmpi(ip_type,'de')&&strcmpi(op_req,'bi')
                op = de2bi(ip,22,'left-msb');
            elseif strcmpi(ip_type,'bi')&&strcmpi(op_req,'de')
                op = bi2de(ip,'left-msb');
            elseif strcmpi(ip_type,'bi')&&strcmpi(op_req,'st')
                tempX = bi2de(ip,'left-msb');
                op = ['x' sprintf('%07d',tempX) 'x'];
            elseif strcmpi(ip_type,'de')&&strcmpi(op_req,'st')
                op = ['x' sprintf('%07d',ip) 'x'];
            end
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
        
        function [nll_app,aic_app,aicc_app,bic_app] = ddm_cost_pdf_nll(obj,x,p)
            if not(isempty(x))
                p = obj.px2p(obj.s.xl,p,x);
            end
            p_RT_and_accuracy = nan(height(obj.data),1);
            
            %while I like this, it takes 50ms
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
            
            for ix_p_config = 1:height(p_mat_unique)
                px = table2struct(p_mat_unique(ix_p_config,:));
                px_array = table2array(p_mat_unique(ix_p_config,:));
                
                [pdf_cw,pdf_cr,~,cdf_cw,cdf_cr] = obj.ddm_pdf(px,obj.s.dt,obj.s.T,obj.s.ddx);
                p_cr = cdf_cr(end);
                p_cw = cdf_cw(end);
                
                %accuracy coding at the moment... should make this flexible
                
                case_config = all(p_mat_array==px_array,2);
                
                
                ix_cr = round(obj.data.rt(case_right&case_config&not(case_nan))/obj.s.dt);
                ix_cw = round(obj.data.rt(case_wrong&case_config&not(case_nan))/obj.s.dt);
                
                pRT_g_cr = pdf_cr(ix_cr)';
                pRT_g_cw = pdf_cw(ix_cw)';
                
                p_RT_and_accuracy(case_right&case_config&not(case_nan)) = pRT_g_cr*p_cr;
                p_RT_and_accuracy(case_wrong&case_config&not(case_nan)) = pRT_g_cw*p_cw;
            end
            
            p_RT_and_accuracy(isnan(p_RT_and_accuracy)) = [];
            p_RT_and_accuracy(p_RT_and_accuracy == 0) = 1e-32;%not great
            
            ll_app = sum(log(p_RT_and_accuracy));
            %     if isnan(ll_app)||isinf(ll_app),error('non scallr ll');end
            if isnan(ll_app),ll_app=-inf;end
            k = length(obj.s.fit_n)+1;%number of free params + 1 for noise
            n = length(p_RT_and_accuracy);
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
            
            p_ = 'a';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 0.6;g_std = 0.15;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_std);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 1.0;
            plbound_.(p_) = 0.01;
            pubound_.(p_) = 7.5;
            prior_.(p_) = @(x) gampdf(x,A_shape,B_scale);
            
            p_ = 't';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 0.3;g_std = 0.075;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_std);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 0.25;
            plbound_.(p_) = 0.1;
            pubound_.(p_) = 0.75;
            prior_.(p_) = @(x) gampdf(x,A_shape,B_scale);
            
            p_ = 'v';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 3;g_std = 1;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_std);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 0.0;
            plbound_.(p_) = -7.5;
            pubound_.(p_) = 7.5;
            prior_.(p_) = @(x) gampdf(x,A_shape,B_scale);
            
            p_ = 'st';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_lo = 0;
            g_up = 0.15;
            pran_.(p_) = unifrnd(g_lo,g_up);
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 0.5;
            prior_.(p_) = @(x) unifpdf(x,g_lo,g_up);
            
            p_ = 'sx';
            modelkey_var{ix} = (p_);ix = ix+1;
            pran_.(p_) = 0;            %not implemented
            pdef_.(p_) = 0;
            plbound_.(p_) = 0;
            pubound_.(p_) = pubound_.a/4;%n.b. the hard reference to a!
            prior_.(p_) = @(x) unifpdf(x,0,pubound_.(p_));
            
        end
    end
    methods (Static)
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_nv(p,dt,T)
            linspace_t = 0:dt:T-dt;
            rt = linspace_t(1:end-1)+dt/2;
            warning('need to think about if this is the right covnersion');
            z = p.x0 + p.a;
            a = 2*p.a;
            v = p.v;
            
            err = 1e-3;%not sure what this should actually be
            tt = rt/(a^2);
            w = z/a;
            p = arrayfun(@(tt) ftt_01w(tt,w,err),tt);
            
            p = p.*exp(-v*a*w - (v^2)*rt/2)/(a^2);
            p_ups = (exp(-2*a*z*v) - 1) / (exp(-2*a*v) - 1);
            %not sure how this is dealing with error trials
            
            pdf_ups = p*p_ups;
            pdf_dow = (1-p)*p_ups;%think this is gennerally wrong...
            cdf_ups = cumtrapz(rt,pdf_ups);
            cdf_dow = cumtrapz(rt,pdf_dow);
            [cdf_ups,cdf_dow] = cdf_t_st(cdf_ups,cdf_dow,rt,t,st,dt);
            
            pdf_ups = diff(cdf_ups)/dt;
            pdf_dow = diff(cdf_dow)/dt;
            %does not currently deal with sv,, sz
            if not(p.sx==0),error('not dealing with this');end
            if not(p.sv==0),error('not dealing with this');end
        end
        function p = ftt_01w(tt,w,err)
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
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf(p,dt,T,ddx)
            %
            
            % %n.b. a multiplication by p.th could stabilise
            x0 = 0;
            %
            linspace_t = 0:dt:T-dt;
            rt = linspace_t(1:end-1)+dt/2;
            %
            xmax = 1.25*p.a+0.2*p.s + p.s;
            xmin = -xmax;
            dx = (xmax-xmin)/ddx;%n.b. if you change the resolution here (100) - then you need to recompile the mex
            xz = [xmax:-dx:xmin]';
            xvm_probe = repmat(xz',length(xz),1)';
            xvm_prev = repmat(xz',length(xz),1);
            % This  has to be different...
            [~,zeroStateIx] = min(abs(xz-(x0)));
            if p.sx == 0
                z0 = zeros(length(xz),1);
                z0(zeroStateIx,1) = 1;
            else
                error('not sure if should be gaussian or uniform');
                %     	z0 = normpdf(xz,xz(zeroStateIx),p.sx);
                z0 = unifpdf(xz,xz(zeroStateIx)-p.s_x0,xz(zeroStateIx)+p.sx);
                z0 = z0/sum(z0);
            end
            %
            sig = p.s*sqrt(dt);
            N_t = length(linspace_t);
            % CORE
            pMat = zeros(length(xz),N_t);
            %
            zn = z0*p.v;
            pMat(:,1) = zn;
            %
            for ix_t = 2:N_t
                xvm_expect = xvm_prev + p.v*dt;
                A = (1/sqrt(2*pi*(sig^2)))*exp(-((xvm_probe - xvm_expect).^2)/(2*(sig^2)));
                An = A./repmat(sum(A,1),length(xz),1);
                An(isnan(An))=0;
                
                An(:,xz>p.a) = 0;
                e = eye(sum(xz>p.a));
                An(1:length(e),1:length(e)) = e;
                An(:,xz<-p.a) = 0;
                e = eye(sum(xz<-p.a));
                An(end-length(e)+1:end,end-length(e)+1:end) = e;
                
                zn = An*zn;
                pMat(:,ix_t) = zn;
            end
            
            cdf_ups = sum(pMat(xz>p.a,:));
            cdf_dow = sum(pMat(xz<-p.a,:));
            %
            [cdf_ups,cdf_dow] = cdf_t_st(cdf_ups,cdf_dow,rt,p.t,p.st,dt);
            
            cdf_ups_end = cdf_ups(end);
            cdf_dow_end = cdf_dow(end);
            cdf_ups = cdf_ups/(cdf_ups_end+cdf_dow_end);
            cdf_dow = cdf_dow/(cdf_ups_end+cdf_dow_end);
            %
            pdf_ups = diff(cdf_ups)/dt;
            pdf_dow = diff(cdf_dow)/dt;
        end
        
        function [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_bru(p,dt,T,N_its)
            
            [vec_rt,vec_correct,~,~,linspace_t] = ddm_def.ddm_bru_core(p,dt,T,N_its);

            cdf_ups = ksdensity(vec_rt(vec_correct),linspace_t,'Support','Positive','function','cdf');
            cdf_dow = ksdensity(vec_rt(not(vec_correct)),linspace_t,'Support','Positive','function','cdf');
            p_ups = sum(vec_correct)/length(vec_correct);
            p_dow = 1-p_ups;
            %this is already handled at the moment.
            
            cdf_ups = cdf_ups*p_ups;
            cdf_dow = cdf_dow*p_dow;
            %
            pdf_ups = diff(cdf_ups)/dt;
            pdf_dow = diff(cdf_dow)/dt;
            rt = linspace_t(1:end-1)+dt/2;
        end
        
        function [vec_rt,vec_correct,x,td_tot,linspace_t] = ddm_bru_core(p,dt,T,N_its)
            
            %%
            linspace_t = 0:dt:T-dt;
            rt = linspace_t(1:end-1)+dt/2;
            %%
            maxIterations = floor(T/dt);
            if length(linspace_t)~=maxIterations,error('Messed up time code');end
            %%
            x0 = 0;
            x0_trial = 2*(rand(N_its,1)-0.5)*(p.sx) + x0;%n.b. this one is uniform.
            x_noise = randn(N_its,maxIterations)*(p.s)*sqrt(dt);
            x_drift = repmat(p.v,N_its,maxIterations)*dt;
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
            
            x_threshold = -p.a;
            [~,x_lower]=sort(x<x_threshold,2,'descend');% should the sign here be the same????
            x_lower=x_lower(:,1);
            x_lower(x_lower==1) = nan;
            
            x_bounds = [x_upper,x_lower];
            x_bounds = gather(x_bounds);
            [ix_time,vec_correct] = nanmin(x_bounds,[],2);
            ix_time(isnan(ix_time)) = length(linspace_t);%not sure if this is best way to deal with problem...
            vec_rt = linspace_t(ix_time)';
            % don't think the time shift can change winner/loser so do it here.
            
            td_tot = p.t + 2*(rand(N_its,1)-0.5)*p.st;
            vec_rt = vec_rt + td_tot;
            
            vec_correct = vec_correct==1;
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
    end
    
    
end

function [cdf_ups,cdf_dow] = cdf_t_st(cdf_ups,cdf_dow,rt,t,st,dt)
td0_vec = find((rt>(t)-st)&(rt<(t)+st));
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



