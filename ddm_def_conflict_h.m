classdef ddm_def_conflict_h < ddm_def
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ddm_def_hist(obj)
            %ovewrite model class property
            obj.modelclass = 'conflict_h';
        end
        
        function get_data(obj)
            get_data@ddm_def(obj);
            %add history for this subject
            obj.data.stim_conflict_hist = [false;obj.data.stim_conflict(1:end-1)];
        end
        
        function p_mat = ddm_cost_add_stim_dependencies(obj,p_mat)
            p_mat.c = obj.data.stim_conflict;
            p_mat.c_hist = obj.data.stim_conflict_hist;
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance(obj)
            
            %use the base method first (rather than redefining the whole thing)
            [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance@ddm_def(obj);
            ix = length(modelkey_var)+1;
            
            p_ = 'bhc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 1;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 6;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'bhnc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 1;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 6;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'vhc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 0.5;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 3;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'vhnc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 0.5;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 3;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'xbc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 0.5;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 3;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            p_ = 'xbnc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 0.5;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 3;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
        end
    end
    
    methods (Static)
        
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
                xvm_expect = xvm_prev + p.v*(...
                    1+p.c*(p.b + p.c_hist*p.bhc + (1-p.c_hist)*p.bhnc)*t_math(ix_t-1)...
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
        
    end
end