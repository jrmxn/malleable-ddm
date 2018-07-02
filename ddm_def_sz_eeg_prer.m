classdef ddm_def_sz_eeg_prer < ddm_def_sz_eeg
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ddm_def_sz_eeg_prer(obj)
            %ovewrite model class property
            obj.modelclass = 'sz_eeg_prer';
            obj.path_data = fullfile('testing','testing_sz.csv');
            obj.info.difficulties = [-5:5];
            obj.info.name_channel = {'nlrt','prer_beta_P3','prer_SSVEP_Cz','prer_beta_Cz'};
        end
        
        
    end
    
    methods (Access = protected)
        
    end
    
    methods (Static)
        
    end
end