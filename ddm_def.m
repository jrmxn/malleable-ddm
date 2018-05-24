classdef ddm_def
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        subject = '';
        modelclass = '';
        id_model = -1+2^4;
        id_fit = 1;
        s = [];%settings
%         init = [];
        modelKey = [];
        %         xl = [];
        data = [];
        opt = [];
        p = [];
%         p_bound = [];
        path_data = '';
    end
    
    methods
        function obj = ddm_def(modelclass)
            obj.modelclass = modelclass;
            obj.modelKey = get_modeldef(obj, 'keyf');
        end
        
        function opt = opt_init(obj)
            opt.TolX = 1e-5;
            opt.MaxIter = 1000;
            opt.parallelsearch = true;
            opt.ps_AccelerateMesh = true;%should only do if smooth
            opt.computeAlgo = 'PS';
        end
        
        function obj = ddm_search(obj)
            if strcmpi(obj.opt.computeAlgo,'PS')
                obj.opt.minoptions = optimoptions(@patternsearch,'Display','iter',...
                    'MaxIter',obj.opt.MaxIter,'TolFun',obj.opt.TolX,'TolX',obj.opt.TolX,...
                    'MaxFunEvals',obj.opt.MaxIter*2,'StepTolerance',1e-3,...
                    'InitialMeshSize',2,'AccelerateMesh',obj.opt.ps_AccelerateMesh,...
                    'UseParallel',obj.opt.parallelsearch,'UseCompletePoll',obj.opt.parallelsearch);%,
            end
            
            if strcmpi(obj.s.minAlgo,'nll')
                ddm_cost_func = @ddm_cost_pdf_nll;
            else
                error('minAlgo not defined');
            end
            
            
            
            obj.opt.init.nll = ddm_cost_func([],obj.opt.init.p,obj.data,obj.s);
            x = p2x(obj.s.xl,obj.opt.init.p);
            x_lb = p2x(obj.s.xl,obj.opt.init.p_lb);
            x_ub = p2x(obj.s.xl,obj.opt.init.p_ub);
            
            if strcmpi(obj.opt.computeAlgo,'PS')
                [x,Fps] = patternsearch(@(x)...
                    ddm_cost_func(x,obj.opt.init.p,obj.data,obj.s)...
                    ,x,[],[],[],[],x_lb,x_ub,obj.opt.minoptions);
            else
                error('Invalid compute algo')
            end
            1;
            obj.p = px2p(obj.s.xl,obj.opt.init.p,x);
        end
        

        function obj = ddm_init(obj, id_model,id_fit)

            obj.s.minAlgo = 'nll';
            obj.s.reinit = false;
            obj.s.dt = 1e-3;
            obj.s.ddx = 100;
            obj.s.T = 5;
            obj.s.inittype = 'random';
            obj.s.path_data = '';
            
            obj.id_model = id_model;
            obj.id_fit = id_fit;
            
            clear v id_model id_fit modelclass subject;
            
            init_p_full = obj.get_modeldef(['init_' obj.s.inittype]);
            
            obj.s.fit_n = sum(obj.debi_model(obj.id_fit,'de','bi'));
            
            
            id_model_index = find(obj.debi_model(obj.id_model,'de','bi'));
            id_fit_index  = find(obj.debi_model(obj.id_fit,'de','bi'));
            
            % Set the parameters over which to optimise from modelType spec
            c = 1;
            for ix_id_fit_index = 1:length(id_fit_index)
                parameter_string = obj.modelKey{id_fit_index(ix_id_fit_index)};
                obj.s.xl.(parameter_string) = c;c = c+1;
            end
            
            % Now set p, using px values (which are themselves either default or loaded)
            for ix_parameter_cell = 1:length(obj.modelKey)
                parameter_string = obj.modelKey{ix_parameter_cell};
                if any(strcmp(parameter_string,obj.modelKey(id_model_index)))
                    init_p_reduced.(parameter_string) = init_p_full.(parameter_string);
                else
                    init_p_reduced.(parameter_string) = 0;
                end
            end
            
            if not(isempty(obj.subject))
                data_ = readtable(obj.path_data);
                %                 data.subject = cell2mat(data.subject);
                data_ = data_(strcmpi(data_.subject,obj.subject),:);
                if not(height(data_)>0),error('Bad subject spec?');end
                obj.data = data_;
                
            end
            
            
            obj.opt = obj.opt_init;
            obj.opt.init.p = init_p_reduced;
            obj.opt.init.p_lb = obj.get_modeldef('lbound');
            obj.opt.init.p_ub = obj.get_modeldef('ubound');
        end
        
        %         function obj = set_subject(obj, subject)
        %             obj.sub = subject;
        %         end
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
        
        
        function outputArg = get_modeldef(obj, deftype)
            ix = 1;
            
            modelkey_var{ix} = 's';ix = ix+1;
            pran_.s = 1;
            pdef_.s = 1;
            plbound_.s = 1;
            pubound_.s = 1;
            
            modelkey_var{ix} = 'a';ix = ix+1;
            g_mea = 0.6;g_std = 0.15;
            [A_shape,B_scale] = gamma_convert(g_mea,g_std);
            pran_.a = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.a = 1.0;
            plbound_.a = 0.01;
            pubound_.a = 7.5;
            
            modelkey_var{ix} = 't';ix = ix+1;
            g_mea = 0.3;g_std = 0.075;
            [A_shape,B_scale] = gamma_convert(g_mea,g_std);
            pran_.t = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.t = 0.25;
            plbound_.t = 0.1;
            pubound_.t = 0.75;
            
            modelkey_var{ix} = 'v';ix = ix+1;
            g_mea = 3;g_std = 1;
            [A_shape,B_scale] = gamma_convert(g_mea,g_std);
            pran_.v = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.v = 0.0;
            plbound_.v = -7.5;
            pubound_.v = 7.5;
            
            modelkey_var{ix} = 'b';ix = ix+1;
            g_sd = 2.5;
            pran_.b = abs(normrnd(0,g_sd));
            pdef_.b = 0.0;
            plbound_.b = 0;
            pubound_.b = 20;
            
            modelkey_var{ix} = 'xb';ix = ix+1;
            g_mea = 0.1;g_std = 0.06;
            [A_shape,B_scale] = gamma_convert(g_mea,g_std);
            pran_.xb = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.xb = 0.0;
            plbound_.xb = 0;
            pubound_.xb = [20];%could def be changed
            
            modelkey_var{ix} = 'st';ix = ix+1;
            g_lo = 0;
            g_up = 0.15;
            pran_.st = unifrnd(g_lo,g_up);
            pdef_.st = 0.0;
            plbound_.st = 0;
            pubound_.st = 0.5;
            
            modelkey_var{ix} = 'sx';ix = ix+1;
            pran_.sx = 0;            %not implemented
            pdef_.sx = 0;
            plbound_.st = 0;
            pubound_.st = plbound_.a(end)/2;
            
            ix = ix;clear ix;
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
            end
        end
        
        
    end
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

function [nll_app,aic_app,aicc_app,bic_app] = ddm_cost_pdf_nll(x,p,data,s)
if not(isempty(x))
    p = px2p(s.xl,p,x);
end
p_RT_and_accuracy = nan(height(data),1);

%while I like this, it takes 50ms
p_mat = struct2table(repmat(p,height(data),1));
p_mat.c = data.stim_conflict;

p_mat_unique = unique(p_mat);
p_mat_array = table2array(p_mat);

for ix_p_config = 1:height(p_mat_unique)
    px = table2struct(p_mat_unique(ix_p_config,:));
    px_array = table2array(p_mat_unique(ix_p_config,:));
    
    [pdf_cw,pdf_cr,t,cdf_cw,cdf_cr] = ddm_pdf(px,s.dt,s.T,s.ddx);
    p_cr = cdf_cr(end);
    p_cw = cdf_cw(end);
    
    
    %accuracy coding at the moment... should make this flexible
    case_right = logical(data.choice);
    case_wrong = not(case_right);
    case_config = all(p_mat_array==px_array,2);
    case_nnan = not(isnan(data.rt));
    %     case_config =
    % %     B = rowfun(@hypot,A,'OutputVariableNames','z')
    
    %     case_incorrect = data.sub_choice=='incorrect';
    %     case_cong_incong = dataSummary.StimVar==px.ip;
    
    if size(t,1)>size(t,2)
        error('vec_t has to be a row vector');
    end
    
    
    %     X = dataSummary.(case_quest).extra.full.(caseCond_ca){1}-t;
    %     [~,ix_ca] = min((abs(X)),[],2);
    ix_cr = round(data.rt(case_right&case_config&case_nnan)/s.dt);
    ix_cw = round(data.rt(case_wrong&case_config)/s.dt);
    
    pRT_g_cr = pdf_cr(ix_cr)';
    pRT_g_cw = pdf_cw(ix_cw)';
    
    p_RT_and_accuracy(case_right&case_config&case_nnan) = pRT_g_cr*p_cr;
    p_RT_and_accuracy(case_wrong&case_config&case_nnan) = pRT_g_cw*p_cw;
end

p_RT_and_accuracy(isnan(p_RT_and_accuracy)) = [];
p_RT_and_accuracy(p_RT_and_accuracy == 0) = 1e-32;%not great

ll_app = sum(log(p_RT_and_accuracy));
%     if isnan(ll_app)||isinf(ll_app),error('non scallr ll');end
if isnan(ll_app),ll_app=-inf;end
k = length(s.fit_n)+1;%number of free params + 1 for noise
n = length(p_RT_and_accuracy);
%prob an issue here with thr fact that some trials are kicked out..
bic_app = log(n)*k-2*ll_app;
aic_app = 2*k-2*ll_app;
aicc_app = aic_app + (2*(k^2) + 2*k)/(n-k-1);
nll_app = -ll_app;
end


function  [pdf_dow,pdf_ups,t_math,cdf_dow,cdf_ups] = ddm_pdf(p,dt,T,ddx)
%

% %n.b. a multiplication by p.th could stabilise
x0 = (-(2*p.c-1)*p.xb);
%
linspace_t = 0:dt:T-dt;
t_math = linspace_t(1:end-1)+dt/2;
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
    %             else
    % 	z0 = normpdf(xz,xz(zeroStateIx),p.sx);
    %                 z0 = unifpdf(xz,xz(zeroStateIx)-p.s_x0,xz(zeroStateIx)+p.sx);
    %                 z0 = z0/sum(z0);
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
    xvm_expect = xvm_prev + p.v*(...
        1+p.c*(p.b)*t_math(ix_t-1)...
        )*dt;
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
td0_vec = find((t_math>(p.t)-p.st)&(t_math<(p.t)+p.st));
if length(td0_vec)<=1
    ix_t0_shift = round((p.t)/dt);
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
%
%             if p.lapser~=0
%                 lapser_slope = linspace_t*p.lapser;
%                 cdf_ups = cdf_ups + lapser_slope;
%                 cdf_dow = cdf_dow + lapser_slope;
%             end
%
cdf_ups_end = cdf_ups(end);
cdf_dow_end = cdf_dow(end);
cdf_ups = cdf_ups/(cdf_ups_end+cdf_dow_end);
cdf_dow = cdf_dow/(cdf_ups_end+cdf_dow_end);
%
pdf_ups = diff(cdf_ups)/dt;
pdf_dow = diff(cdf_dow)/dt;
end

function [A_shape,B_scale] = gamma_convert(g_mea,g_std)
alpha_shape = (g_mea^2)/(g_std^2);
beta_rate = (g_mea)/(g_std^2);
A_shape = alpha_shape;
B_scale = 1/beta_rate;
end