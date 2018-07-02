classdef ddm_def_sz_eeg_pos < ddm_def_sz_eeg
    %DDM_DEF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ddm_def_sz_eeg_pos(obj)
            %ovewrite model class property
            obj.modelclass = 'sz_eeg_pos';
            obj.path_data = fullfile('testing','testing_sz.csv');
            obj.info.difficulties = [-5:5];
            obj.info.name_channel = {'nlrt','pos_random_O','pos_alpha_O','pos_ssvep_O','pos_theta_O'};%,'prer_beta_P3','prer_SSVEP_Cz','prer_SSVEP_Cz'};
        end
        
        
    end
    
    methods (Access = protected)
        
    end
    
    methods (Static)
        
    end
end