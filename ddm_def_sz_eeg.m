classdef ddm_def_sz_eeg < ddm_def_sz
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ddm_def_sz_eeg(obj)
            %ovewrite model class property
            obj.modelclass = 'sz_eeg';
            obj.path_data = fullfile('testing','testing_sz.csv');
            obj.info.difficulties = [-5:5];
        end
        
        function get_data(obj)
            get_data@ddm_def_sz(obj);
        end
        
        function p_mat = ddm_cost_add_stim_dependencies(obj,p_mat)
            p_mat = ddm_cost_add_stim_dependencies@ddm_def(obj,p_mat);
            p_mat.difficulty = obj.data.difficulty;
            p_mat.pre_alpha_O = obj.data.pre_alpha_O;
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance(obj)
            function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = def_eeg_params(p_, g_sd)
                modelkey_var = (p_);
                pd_hn_i = makedist('Normal','mu',0,'sigma',g_sd);
                pran_ = pd_hn_i.random;
                pdef_ = 0.0;
                plbound_ = -5*g_sd;
                pubound_ = +5*g_sd;
                prior_ = @(x) pdf(pd_hn_i,x);
            end
            
            [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ddm_def_instance@ddm_def_sz(obj);
            ix = length(modelkey_var)+1;
            
            % we can put eeg based noise on...
            p_ = 't_pre_alpha_O';
            [modelkey_var{ix},pran_.(p_),pdef_.(p_),plbound_.(p_),pubound_.(p_),prior_.(p_)] ...
                = def_eeg_params(p_, 0.1);ix = ix+1;
            p_ = 'z_pre_alpha_O';
            [modelkey_var{ix},pran_.(p_),pdef_.(p_),plbound_.(p_),pubound_.(p_),prior_.(p_)] ...
                = def_eeg_params(p_, 0.1);ix = ix+1;
            p_ = 'v_pre_alpha_O';
            [modelkey_var{ix},pran_.(p_),pdef_.(p_),plbound_.(p_),pubound_.(p_),prior_.(p_)] ...
                = def_eeg_params(p_, 0.1);ix = ix+1;
        end
    end
    
    methods (Static)
        
        function  pdf_ = ddm_prt_ana(p, rt, eeg_mod)
            err = 1e-8;
            
            diffi_str = ddm_def_sz.diff2drift(p.difficulty);
            px = p;
            %n.b. the re-assignment to v here is also critical for v_eeg interactions
            px.v = px.(diffi_str);
            for ix_eeg_mod = 1:length(eeg_mod)
                ch_str = eeg_mod(ix_eeg_mod).channel;
                p_str = eeg_mod(ix_eeg_mod).param;
                px.(p_str) = px.(p_str) + ...
                    px.(ch_str) * px.(sprintf('%s_%s',p_str,ch_str));
            end
            
            h_pdf = @(x) ddm_def.hddm_pdf_full(x,px.v,px.sv,px.a,px.z,px.sz,px.t,px.st,err);
            pdf_ = arrayfun(@(x) h_pdf(x),+rt);
            
        end
        
    end
end