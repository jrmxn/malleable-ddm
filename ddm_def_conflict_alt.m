classdef ddm_def_conflict_alt < ddm_def
    
    properties
    end
    
    methods
        function obj = ddm_def_conflict_alt
            %light initialisation so functions can be used easily
            obj.modelclass = 'conflict_alt';
            obj.path_data = fullfile('testing','testing.csv');
        end
        
        function p_mat = ddm_cost_add_stim_dependencies(obj,p_mat)
            p_mat.c = obj.data.stim_conflict;
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance(obj)
            
            %use the base method first (rather than redefining the whole thing)
            [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance@ddm_def(obj);
            ix = length(modelkey_var)+1;
            
            p_ = 'zc';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_alpha = 1;
            g_beta = 3;
            pd_hn = makedist('beta','a',g_alpha,'b',g_beta);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = -1;
            pubound_.(p_) = 1;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            for ix_b = 1:9
            p_ = sprintf('b%d',ix_b);
            modelkey_var{ix} = (p_);ix = ix+1;
            g_sd = 5;
            pd_hn = makedist('HalfNormal','mu',0,'sigma',g_sd);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 20;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            end
        end
        
        
    end
    methods (Static)
        
        function [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_bru(p,lt,N_its)
            
            %             if not(p.sz==0)
            %                 error('Model not defined for sz not equal to zero');
            %             end
            rt = (lt(1:end-1)+lt(2:end))*0.5;
            
            dt = lt(2)-lt(1);
            T = lt(end)+dt;
            maxIterations = floor(T/dt);
            if (length(lt)-1)~=maxIterations,error('Messed up time code');end
            %%
            %modification for conflict
            z = p.z - 0.5*(2*p.c-1)*p.zc;
            za = z*p.a;
            za_trial = (rand(N_its,1)-0.5)*(p.sz) + za;%n.b. this one is uniform.
            x_noise = randn(N_its,maxIterations)*(p.s)*sqrt(dt);
            
            %modification for conflict
            V = repmat((p.v + p.sv*randn(N_its,1)),1,maxIterations);
            CB = repmat(p.c*p.b1*rt,N_its,1);
            %             x_drift = V
            %%
            x = nan(N_its,maxIterations);
            x(:,1) = za_trial;
            for ix_t = 2:maxIterations
                error('not implemented');
%                 y2 = -sign(z-p.z)*(x(:,ix_t-1)-za_trial);
%                 y3 = y2./(1+exp(-25*y2));
%                 y4 = -sign(z-p.z)*(x(:,ix_t-1)-p.z*p.a);
%                 y5 = y4./(1+exp(-25*y4));
                dx = x_noise(:,ix_t) + V(:,ix_t).*( ...
                    +1 ...
                    +CB(:,ix_t) ...
                    +p.b2*y2 ...
                    +p.b3*y3 ...
                    +p.b4*y4 ...
                    +p.b5*y5 ...
                    )*dt;
                x(:,ix_t) = x(:,ix_t-1) + dx;
                %
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
        end
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_trm(p,lt,dx)
            %
            
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
            z = p.z - 0.5*(2*p.c-1)*p.zc;
            za = (z)*p.a;
            
            [~,zeroStateIx] = min(abs(xz-za));
            if (p.sz <= 1e-3)
                x0 = zeros(length(xz),1);
                x0(zeroStateIx,1) = 1;
            else
%                 error('Model not defined for sz not equal to zero');
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
                y2 = -sign(z-p.z).*(xz-za);
                y3 = y2./(1+exp(-100*y2));
                y4 = -sign(z-p.z).*(xz-p.z*p.a);
                y5 = y4./(1+exp(-100*y4));
                    xvm_expect = xvm_prev + v*(...
                        1 ...
                        +p.b1*p.c*lt(ix_t-1)...
                        +p.b2*y2...
                        +p.b3*y3...
                        +p.b4*y4...
                        +p.b5*y5...
                        )*dt;
                    
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




