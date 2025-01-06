%% this code investigates the effect of varrying network parameters
% one at a time, PN-KC weights, KCs spiking thresholds and number of inputs
% on the memory performance of model fruit-flies Drosophila
% this simulation uses realistic odor responses recordings from
% (Hallem-Carlson 2006) and transformed from ORN to PN responses as by
% Olsen et.al; detailed in the Methods section in our paper.
%%

clear all;
clc
% param. initializations
numModelInstances = 20; % trials, each corresponding to different mushroom body instance

numTrials = 30; %number of noisy trials per each odors

%loading Hallem-Carlson data
[PNactivity,~,~] = propagateORNs2PNs;
PNactivity = PNactivity';
K_odors = size(PNactivity,2);
nPNs = size(PNactivity,1);

% levels of noise added to inputs
% NScales=3;
noiseLevels = [0.5,1,2];

MBmodels_allVar_noHomeo = {};
MBmodels_allVar_withHomeo = {};
allTestedOdors = {};

% the odors we used in Fig2, 3 and 4, stored in this .mat file.
% load('Data_submitted_fly_wNoise11.mat','PNtrials');

%% loops over different number of input odors, N=50,100,150,200.
% testing models at different number of odors, Fig2C1

% random fly instantiations
for randomTrials=1:numModelInstances
    
    % MBmodel.nKCs = 2000; % number of neurons in the hidden layer/ KCs
    % MBmodel.nPNs=24;  %number of dimensions in the input data/ PNs

    % new model instantiation function
    % buildMBmodel(nKCs, nPNs ,isVar_Nclaws, isVar_weights, isVar_Theta, varargin)
    % ->varargin for mean, std and lower/upper bounds of Nclaws and theta
    MBmodel_varN_varW_varTh = MBmodelBuilder(2000,24, true, true, true); %formerly {name} without special extension
    % MBmodel_fixN_fixW_varTh = MBmodelBuilder(2000,24, false,false,true); % formerly {name}_Homog
    % MBmodel_varN_varW_fixTh = MBmodelBuilder(2000,24, true, true, false, 'theta_mean',10); % formerly {name}_fixedTheta
    % MBmodel_fixN_fixW_fixTh = MBmodelBuilder(2000,24, false,false,false, 'theta_mean',10+randn(1)); % formerly {name}_Homog_Ftheta
    % ... etc 

    noiseScale = 2;
    noise = noiseLevels(noiseScale)  % noiseLevels=[0.5,1,2];
    
    PNtrials = zeros(nPNs, K_odors, numTrials);
    PNtrials(:,:,1) = PNactivity; %first trial is "simply" the first instance of the fictitious odor response
    % for all trials, add some noise to the generated fictitious odor response, and calculate a SNR ratio
	for t = 1:numTrials-1
		% actually, why not sample the randn once and use the SAME value for PNtrials and SNR?
        PNtrials(:,:,t+1) = PNactivity + ... % add some noise to the generated fictitious odor response
            getPNStdevBhandawat(PNactivity) .* ... %returns hard-coded glomeruli-specific std-value from Bhandawat2007
            noise.*randn(nPNs, K_odors); %noise amplitude times normal distr. number
        % SNR(:,:,t)= (PNactivity./(getPNStdevBhandawat(PNactivity) .* ...
        %     noise.*randn(nPNs, K_odors))).^2; % an individual SNR ratio (x is signal, glomeruli-spec std times random noise amplitude
        % SNR(:,:,t)=log(SNR(:,:,t)); %log-SNR 
    end
        
    % making sure that PN responses are +ve
    PNtrials(PNtrials<0)=0;
    
    %rescaling the PN responses in range from 0 to 5
	%  --> WHY?
    PNtrials = rescale(PNtrials,0,5);
    PNtrials = reshape(PNtrials, size(PNtrials,1),[]); %make 2-D
    
    %start the parallel pool
	%% begin optimisation procedure to get C_theta and APLgain values that satisfy coding level constraints
	% for the variable weight, variable theta model.

    %optimises and returns an updated MBmodel
	MBmodel_varN_varW_varTh_noHomeo = optimiseMBparams_APLgain_Ctheta(MBmodel_varN_varW_varTh, PNtrials);
    MBmodel_varN_varW_varTh_withHomeo = optimiseMBparams_homeostaticExcitation(MBmodel_varN_varW_varTh, PNtrials);
    % theta = (MBmodel_varN_varW_varTh.C_theta.*MBmodel_varN_varW_varTh.theta) ;
 
    MBmodel = MBmodel_varN_varW_varTh_withHomeo;
    MBmodel.APLgain = MBmodel.alpha; % just make an alias

    %% save the parameters for each fly, connectivity weights, spiking thresholds, inhibition Gains,for all 8 models
    
    % save the flies tuned on original Hallem-Carlson Data:        
    MBmodels_allVar_noHomeo{end+1} = MBmodel_varN_varW_varTh_noHomeo;
    MBmodels_allVar_withHomeo{end+1} = MBmodel_varN_varW_varTh_withHomeo;
    allTestedOdors{end+1} = PNtrials;
			
end

% 
% %% calculate KC responses (Y) to all odors
% 
% Y = calculateKCresponse(MBmodel,PNtrials);
% totalKCactivity = 
