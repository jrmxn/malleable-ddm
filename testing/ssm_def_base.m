classdef ssm_def_base < ssm_def
    
    properties
    end
    
    methods
        function obj = ssm_def_base
            %light initialisation so functions can be used easily
            obj.modelclass = 'base';
            obj.path_data = fullfile('testing.csv');
        end
        
        function p_mat = ssm_cost_add_stim_dependencies(obj,p_mat)
            % this function could be omitted for the basic ssm fits
            % but adding to show how to incorporate trial conditions
            % (conflict in this case)
            p_mat = ssm_cost_add_stim_dependencies@ssm_def(obj,p_mat);
            p_mat.c = obj.data.stim_conflict;
        end
        
    end
    
    methods (Access = protected)
        function [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ssm_def_instance(obj)
            %use the base method first (rather than redefining the whole thing)
            [modelkey_var,pran_,pdef_,plbound_,pubound_,prior_] = ssm_def_instance@ssm_def(obj);
            ix = length(modelkey_var)+1;
        end
    end
    methods (Static)
        
        
    end
end




