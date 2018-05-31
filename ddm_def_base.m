classdef ddm_def_base < ddm_def
    
    properties
    end
    
    methods
        function obj = ddm_def_base
            %light initialisation so functions can be used easily
            obj.modelclass = 'base';
            obj.path_data = 'testing.csv';
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
        end
    end
    methods (Static)
        
        
    end
end




