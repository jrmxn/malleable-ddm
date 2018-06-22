classdef ddm_def_sz_ext < ddm_def_sz
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ddm_def_sz_ext(obj)
            %ovewrite model class property
            obj.modelclass = 'sz_ext';
            obj.path_data = fullfile('testing','testing_sz.csv');
            obj.info.difficulties = [-5:5];
            %             obj.info.name_channel = {'pre_alpha_O','pre_ssvep_O','pre_random_O','pre_mlfr_C','pre_theta_FC','nlrt'};
        end
        
        function get_data(obj)
            get_data@ddm_def_sz(obj);
        end
        
        function p_mat = ddm_cost_add_stim_dependencies(obj,p_mat)
            p_mat = ddm_cost_add_stim_dependencies@ddm_def_sz(obj,p_mat);
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance(obj)
            
            [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance@ddm_def_sz(obj);
            ix = length(modelkey_var)+1;
            
            p_ = 'k';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 0.5;
            pd_hn = makedist('Normal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = -12;
            pubound_.(p_) = +12;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
        end
    end
    
    methods (Static)
        
        function  [pdf_,p_cr] = ddm_prt_ana(p,rt)
            %this doesn't use leakage - just for initialisation
            if isfield(p,'k')
                if not(p.k==0),error('Non zero leak not allowed with this pdf method');end
            end
            [pdf_, p_cr] = ddm_prt_ana@ddm_def_sz(p,rt);
        end
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups,correct_side] = ddm_pdf_ana(p,rt)
            %this doesn't use leakage - just for initialisation
            if isfield(p,'k')
                if not(p.k==0),error('Non zero leak not allowed with this pdf method');end
            end
            [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups,correct_side] = ddm_pdf_ana@ddm_def_sz(p,rt);
        end
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups,correct_side] = ddm_pdf_trm(p,lt,dx)
            %
            diffi_str = ddm_def_sz.diff2drift(p.difficulty);
            p.v = p.(diffi_str);
            correct_side = sign(p.v);
            
            dt = lt(2)-lt(1);
            f = 5;%obj.s.x_bound_scale;%not sure what the consequence of shrinking this is
            xmax = p.a + f*sqrt(dt)*p.s;
            xmin = 0 - f*sqrt(dt)*p.s;
            
            xz = [xmax:-dx:xmin]';
            xvm_probe = repmat(xz',length(xz),1)';
            xvm_prev = repmat(xz',length(xz),1);
            
            %modification for conflict
            %             za = (p.z*p.a);
            %defined so that p.zc = 1 is the max conflict bias you can have (with z =
            %0.5)
            z = p.z;
            za = (z)*p.a;
            
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
                
                zn = x0;
                xz_ups = xz>p.a;
                xz_dow = xz<0;
                e_ups = eye(sum(xz_ups));
                e_dow = eye(sum(xz_dow));
                len_e_ups = length(e_ups);
                len_e_dow = length(e_dow);
                
                pMat(:,1) = zn;
                for ix_t = 2:N_t
                    %modification for conflict here.
                    y = (p.k)*(xz-za);
                    xvm_expect = xvm_prev + (v-y)*dt;
                    
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
        
    end
end