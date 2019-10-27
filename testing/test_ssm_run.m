% Simplest demonstration of ssm_def
% Define the simplest inheriting function ssm_def_base.m
%%
clear;
f_path_data = 'testing.csv';
addpath('..');
addpath(fullfile('..','ext_mod'));
for ix_fit_type = [1, 2, 3]
    for ix_sub = 1:3
        
        clear sr;
        sr = ssm_def_base;
        sr.subject = sprintf('sub%02d',ix_sub);
        sr.path_data = fullfile('testing.csv');
        mk = sr.ssm_get_instance('keyr');
        rng(ix_sub);
        %% Define the model parameters included in the model
        id_model = sr.debi_model(0,'de','bi');
        id_model(mk.s) = 1;
        id_model(mk.z) = 1;
        id_model(mk.a) = 1;
        id_model(mk.t) = 1;
        id_model(mk.v) = 1;
        id_model(mk.st) = 1;
        id_model(mk.sz) = 1;
%         id_model(mk.sv) = 1;
        id_model_de = sr.debi_model(id_model,'bi','de');
        
        %% Define the model parameters to fit
        id_search = sr.debi_model(0,'de','bi');
        id_search(mk.s) = 0;%noise fixed at s=1, do not optimise
        id_search(mk.z) = 0;%starting point fixed 
        id_search(mk.a) = 1;%fit the threshold
        id_search(mk.t) = 1;%fit the non-decision time
        id_search(mk.v) = 1;%fit the drift
        id_search(mk.st) = 1;%fit the noise in non-decision time
        id_search(mk.sz) = 1;%fit noise in starting point
%         id_search(mk.sv) = 1;%fit noise in drift (slow)
        id_search_de = sr.debi_model(id_search,'bi','de');

        % Initialise the fit
        sr.ssm_init(id_model_de,id_search_de);
        
        %temporary
        sr.opt.parallelsearch = false;
        
        
        % Set the method used to calculate the likelihood:
        % In reality you would choose one of these methods, ssm_pdf_ana
        % (converted from Navaro and Fuss and HDDM) if you have a DDM
        % (including if you have noise non-decision time, noise in drift
        % and noise in starting point).
        % For more complex models you would use ssm_pdf_trm which defines
        % the ssm with a transition matrix. For specific cases (and to
        % check correct specification of transition matrix), you may want
        % to compute with a brute force method (ssm_prt_ana) - not really
        % intended for optimisation.
        if ix_fit_type ==1
            %the HDDM method
            sr.modelclass = 'base_ana_prt';%just a name for the folder
            sr.ssm_pdf = @(a,b) sr.ssm_prt_ana(a,b);
        elseif ix_fit_type==2
            %the transition matrix method
            sr.modelclass = 'base_trm';%just a name for the folder
            %use relaxed values to get in ballpark of solution:
            sr.s.dx = 0.025;%speed things up a bit  (not as accurate, for final optimisation would say <0.005)
            sr.s.dt = 1e-3;%speed things up a bit (not as accurate, for final optimisation would say <5e-4)
            sr.ssm_pdf = @(a,b) sr.ssm_pdf_trm(a,b,sr.s.dx);
        elseif ix_fit_type==3
            %brute force method
            sr.modelclass = 'base_bru';%just a name for the folder
            sr.s.nits = 10e3;%(only 10e3 - not very accurate)
            sr.ssm_pdf = @(a,b) sr.ssm_pdf_bru(a,b,sr.s.nits);
        elseif ix_fit_type == 4
            % buggy - needs fixing
            %the HDDM method used to evaluate the entire PDF (a bit slower,
            %mainly used for plotting)
            sr.modelclass = 'base_ana';%just a name for the folder
            sr.ssm_pdf = @(a,b) sr.ssm_pdf_ana(a,b);
        end
        
        sr.ssm_fit;
        sr.ssm_save;
    end
end
