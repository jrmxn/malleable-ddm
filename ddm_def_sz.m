classdef ddm_def_sz < ddm_def
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ddm_def_sz(obj)
            %ovewrite model class property
            obj.modelclass = 'sz';
            obj.path_data = fullfile('testing','testing_sz.csv');
            obj.info.difficulties = [-5:5];
        end
        
        function [p_mat, sr_full] = aux_gather(obj,f_path,id_model_de,id_search_de,sub_cell,minorfin)
            [p_mat, sr_full] = aux_gather@ddm_def(obj,f_path,id_model_de,id_search_de,sub_cell,minorfin);
            
            issz = nan(height(p_mat),1);
            for ix_sub = 1:height(p_mat)
                if isfield(sr_full(ix_sub).info,'issz')
                    issz(ix_sub) = sr_full(ix_sub).info.issz;
                end
            end
            if all((issz==1)|(issz==0)),issz = logical(issz);end
            p_mat.issz = issz;
            %used to be: (but above deals with missing data)
            %p_mat.issz = arrayfun(@(ix) sr_full(ix).info.issz,[1:height(p_mat)])';
        end
        
        function get_data(obj)
            get_data@ddm_def(obj);
            
            %for now just exclude gender disc.
            case_em = strcmpi(obj.data.cue_question,'EM');
            
            obj.data.rt = obj.data.sub_respTime_lsl;
            obj.data.rt(not(case_em)) = nan;
            
            obj.data.choice = strcmpi( obj.data.sub_choice,'HA')-strcmpi( obj.data.sub_choice,'AN');
            obj.data.choice(obj.data.choice<=0) = 0;%choice is 0 or 1, not -1 or 1
            obj.data.choice(not(case_em)) = nan;
            
            obj.data.difficulty = round(obj.data.cue_coh_signed*10);
            
            %reduce the data (less confusing...)... but annoying for
            %inheritance.
            %            obj.data = obj.data(:,{'choice','rt','difficulty'});
        end
        
        function p_mat = ddm_cost_add_stim_dependencies(obj,p_mat)
            p_mat = ddm_cost_add_stim_dependencies@ddm_def(obj,p_mat);
            p_mat.difficulty = obj.data.difficulty;

            % A bit ugly that this is here, but if specific difficulty is not in the
            % model name, then set it to nan in p_mat so that the cost
            % function does not use it.
            case_skip = false(size(p_mat.difficulty));
            model_def = obj.ddm_print_model(true);
            for ix_diff = 1:length(obj.info.difficulties)
                diff_ = obj.info.difficulties(ix_diff);
                diffi_case = contains(model_def,obj.diff2drift(diff_));
                if not(any(diffi_case))
                    case_skip = case_skip|(p_mat.difficulty == diff_);
                end
            end
            p_mat.skip(case_skip) = 1;  
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance(obj)
            
            %redfine the whole thing
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
            
            p_ = 'a';
            modelkey_var{ix} = (p_);ix = ix+1;
            g_mea = 1;g_sd = 0.25;
            [A_shape,B_scale] = obj.gamma_convert(g_mea,g_sd);
            pran_.(p_) = gamrnd(A_shape,B_scale,[1,1]);
            pdef_.(p_) = 1.0;
            plbound_.(p_) = 0.05;
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
            g_beta = 3;
            pd_hn = makedist('beta','a',g_alpha,'b',g_beta);
            pran_.(p_) = pd_hn.random;
            pdef_.(p_) = 0.0;
            plbound_.(p_) = 0;
            pubound_.(p_) = 1;
            prior_.(p_) = @(x) pdf(pd_hn,x);
            
            
            diffi_val = unique(obj.info.difficulties);
            diffi_str = arrayfun(@(x) ddm_def_sz.diff2drift(x),diffi_val,'UniformOutput',false);
            for ix_v = 1:length(diffi_str)
                p_ = diffi_str{ix_v};
                modelkey_var{ix} = (p_);ix = ix+1;
                g_mea = diffi_val(ix_v)/2;g_sd = 0.5;
                pd_hn = makedist('Normal','mu',g_mea,'sigma',g_sd);
                pran_.(p_) = pd_hn.random;
                pdef_.(p_) = 0.0;
                plbound_.(p_) = -15;
                pubound_.(p_) = +15;
                prior_.(p_) = @(x) pdf(pd_hn,x);
            end
            
        end
    end
    
    methods (Static)
        
        function  [pdf_,p_cr] = ddm_prt_ana(p,rt)
            diffi_str = ddm_def_sz.diff2drift(p.difficulty);
            px = p;
            px.v = p.(diffi_str);
            [pdf_, p_cr] = ddm_prt_ana@ddm_def(px,rt);
        end
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups,correct_side] = ddm_pdf_ana(p,rt)
            diffi_str = ddm_def_sz.diff2drift(p.difficulty);
            px = p;
            px.v = p.(diffi_str);
            [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_ana@ddm_def(px,rt);
            correct_side = sign(px.v);
        end
        
        function  [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups,correct_side] = ddm_pdf_trm(p,rt,dx)
            diffi_str = ddm_def_sz.diff2drift(p.difficulty);
            px = p;
            px.v = p.(diffi_str);
            [pdf_dow,pdf_ups,rt,cdf_dow,cdf_ups] = ddm_pdf_trm@ddm_def(px,rt,dx);
            correct_side = sign(px.v);
        end
        
        function diffi_str = diff2drift(diffi_val)
            diffi_str = sprintf('v+%01d',diffi_val);
            %             diff10_str = arrayfun(@(x) sprintf('v+%01d',x),diff10_val,'UniformOutput',false);
            diffi_str = strrep(diffi_str,'+-','m');
            diffi_str = strrep(diffi_str,'+','p');
            diffi_str = strrep(diffi_str,'vp0','vz');
        end
    end
end