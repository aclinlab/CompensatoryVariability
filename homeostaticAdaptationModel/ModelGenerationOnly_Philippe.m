%% this code investigates the effect of varrying network parameters
% one at a time, PN-KC weights, KCs spiking thresholds and number of inputs
% on the memory performance of model fruit-flies Drosophila
% this simulation uses realistic odor responses recordings from
% (Hallem-Carlson 2006) and transformed from ORN to PN responses as by
% Olsen et.al; detailed in the Methods section in our paper.
%%

clear all;
clc
%loading Hallem-Carlson data
load('hallem_olsen.mat');
% param. initializations
numModelInstances = 10; % trials, each corresponding to different mushroom body instance
numtrainingSamples=15;

% at 4 scales of 200 odors, in increment of 50: 50,100,150,200.
% to test at just one scale (like 100 odors): simply change
% modss=1, and mulOd=100
modss=4; %4
mulOd=200; %200 odors we want to test

odors=mulOd;

numTrials = 30; %number of noisy trials per each odors
PN = hallem_olsen(1:110,:)';
PNs_2_=zeros(24,odors-100); %we want to get to mulOd total odors, and will retrieve 100 from another file
k=0;
binsAll=zeros(24,100);

% Sigmoid_derivative = @(x) exp(0.9.*x)./((1.+exp(0.9.*x)).^2)
% Sigmoid_derivative = @(x) exp(-0.9.*x)./((1.+exp(-0.9.*x)).^2)
% Sigmoid_derivative = @(x) exp(-10.*x)./((1.+exp(-10.*x)).^2)

% C_SoftMax=1e-4;
% levels of noise added to inputs
NScales=3;
% to be explored
noiseLevels=[0.5,1,2];

% learning rates for the KC-MBON weights.
lrs=10;
LRs=10.^[-5 -4 -3 -2.75 -2.5 -2.25 -2 -1 0 1];
% scales for indetermincy in the decision making (Softmax function--> Step
% function) 5 scales of C investigated in Fig2C4
Crange=5;

MBmodels_allVar = {};

%% create artificial odors, n odors
for Pn=1:24
    [prob,bins]=hist(PN(Pn,:),100);
    prob=prob/sum(prob);
    
    binsAll(Pn,:) = bins;
    PNs_2_(Pn,k+1:k+odors-100)=randsample(bins,(odors-100),'true',prob);
    
end

% the odors we used in Fig2, 3 and 4, stored in this .mat file.
load('Data_submitted_fly_wNoise11.mat','PNtrials');

% the scaled PNs responses, the base 100 fictitious odors.
PNs_1_=PNtrials(:,:,1);

%% recover the rescaling factors

% get the maximum bin center for each PN derived from the original
% Hallem-Olsen data
maxRespBinPerPN = max(binsAll,[],2);

% get the maximum response in the rescaled randomly resampled PNs
maxRespPNsRescaledPerPN = max(PNs_1_,[],2);

PNsAboveBestFit = true(24,1);

% In this while loop:
% draw a best fit line comparing the maximum response in the rescaled
% randomly resample PNs to the maximum bin center for each PN in the
% original H-O data. Most PNs will match, but in some PNs, by random chance
% they will not have sampled the top response in 100 odors. These PNs will
% be below the best fit line, while the "matching" PNs will be above. On
% the next iteration of the while loop, redraw the best fit using only the
% PNs that lie above the best fit line from the current iteration
% continue until none of the PNs are above the best fit line (because they
% lie on it almost exactly)

while sum(PNsAboveBestFit)
    sum(PNsAboveBestFit)
    p = polyfit(maxRespBinPerPN(PNsAboveBestFit), maxRespPNsRescaledPerPN(PNsAboveBestFit),1);
    
    % the correct PNs will be above the best fit line because the incorrect PNs
    % are outliers that drag down the line of best fit
    PNsAboveBestFit = maxRespPNsRescaledPerPN > ( maxRespBinPerPN*p(1) + p(2) +0.000001);
    % the 0.000001 is for rounding errors, otherwise you end up in endless
    % loops
end

%% plot the data and the best fit line to confirm it passes through the correct PNs
figure,scatter(maxRespBinPerPN,maxRespPNsRescaledPerPN);
hold on
xrange = [min(binsAll(:)) max(binsAll(:))];
plot(xrange, p(1)*xrange + p(2))

%% recover the original fictitious odors from the rescaled fictitious odors
% why dow e do this both ways? i.e. rescale PN_1 to PN_2 scale and PN_2 to
% PN_1? Should we not map onto the same scale by remapping only 1 set?
PNs_1 = (PNs_1_ - p(2))/p(1);
PNs_2=  (PNs_2_*p(1))+p(2);

% Why do we take PN_2_ instead of PN_2, this looks like a typo
%or maybe on purpose, because both PN_1 and PN_2 are rescaled, which seemed
%odd in the first place?
PNs=[PNs_1 PNs_2_]; % now this is the 24 by 200 matrix (24 ORNs 200 odors)

%% loops over different number of input odors, N=50,100,150,200.
% testing models at different number of odors, Fig2C1

% for mods=1:modss
mods = 2
    odors= mulOd*(mods/modss);
    
    % input fictitious odors sampled from Hallem-Carlson
    % this slicing operation selects only the (rescaled) odours from
    % Hallem-Carlson
    x=PNs(:,1:odors);
    
    % original inputs from Hallem-Carlson
    %    x=PN(:,1:odors);
    
    % random fly instantiations
    for randomTrials=1:numModelInstances
        
        % random split of 'good' and 'bad' odors assignments: class 1 Vs.
        % class2
        classAction1=randsample([1:odors],round(odors/2),'false');
        
        % APLgainP= zeros(1,6);
        % APLgains=zeros(1,9);
        
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
                
        % testing models at different scales of input noise: Fig2C3
        %for noiseScale=1:NScales
         noiseScale = 2   
        noise= noiseLevels(noiseScale);  % noiseLevels=[0.5,1,2];
        
        PNtrials = zeros(24, odors, numTrials);
        SNR=zeros(24,odors,numTrials);
        PNtrials(:,:,1) =x; %first trial is "simply" the first instance of the fictitious odor response
        % for all trials, add some noise to the generated fictitious odor response, and calculate a SNR ratio
		for t = 1:numTrials-1
			% actually, why not sample the randn once and use the SAME value for PNtrials and SNR?
            PNtrials(:,:,t+1) = x + ... % add some noise to the generated fictitious odor response
                getPNStdevBhandawat(x) .* ... %returns hard-coded glomeruli-specific std-value from Bhandawat2007
                noise.*randn(24, odors); %noise amplitude times normal distr. number
            SNR(:,:,t)= (x./(getPNStdevBhandawat(x) .* ...
                noise.*randn(24, odors))).^2; % an individual SNR ratio (x is signal, glomeruli-spec std times random noise amplitude
            SNR(:,:,t)=log(SNR(:,:,t)); %log-SNR 
        end
        
		% mean SNR given odor sample size and given MB instance
        SNR_averageAcrossTrials_andPNs(mods,randomTrials,noiseScale)= mean(SNR(:)); 
        
        % making sure that PN responses are +ve
        PNtrials(PNtrials<0)=0;
        
        %rescaling the PN responses in range from 0 to 5
		%  --> WHY?
        PNtrials=rescale(PNtrials,0,5);
        PNtrials = reshape(PNtrials, size(PNtrials,1),[]); %make 2-D
        
        %start the parallel pool
		%% begin optimisation procedure to get C_theta and APLgain values that satisfy coding level constraints
		% for the variable weight, variable theta model.
        
        % CLevelP=[];
        % INHAbs_CLP=[];

        tic
        %optimises and returns an updated MBmodel
		MBmodel_varN_varW_varTh = optimiseMBparams_APLgain_Ctheta(MBmodel_varN_varW_varTh, PNtrials)
        
        % theta = (MBmodel_varN_varW_varTh.C_theta.*MBmodel_varN_varW_varTh.theta) ;
     
	    toc
        MBmodel = MBmodel_varN_varW_varTh;
        MBmodel.APLgain = MBmodel.alpha; % just make an alias

        %% save the parameters for each fly, connectivity weights, spiking thresholds, inhibition Gains,for all 8 models
        
        % if(odors==100)
            % save(strcat([pwd,'2B_C_VarDegradesPerf_fly_wNoise',num2str(randomTrials),num2str(noiseScale)]) , 'thisW_HomogModel','APLgains',...
                % 'thetaH_Ftheta', 'thisW','thetaS_Ftheta' ,'thisW_varWonly','thisW_varN','theta_varw_fixedN_scalar','theta_varN_fixedW_scalar',...
                % 'thetaH', 'thetaS','theta_varw_fixedN', 'theta_varN_fixedw','PNtrials','classAction1' );
            
        % end
        % save the flies tuned on original Hallem-Carlson Data:
        % if(odors==110)
        % save(strcat(['S1_VarDegradesPerf_HOInp_fly_wNoise',num2str(randomTrials),num2str(noiseScale)]) , 'thisW_HomogModel','APLgains',...
        % 'thetaH_Ftheta', 'thisW','thetaS_Ftheta' ,'thisW_varWonly','thisW_varN','theta_varw_fixedN_scalar','theta_varN_fixedW_scalar',...
        % 'thetaH', 'thetaS','theta_varw_fixedN', 'theta_varN_fixedw','PNtrials','classAction1' );
        %
        % end
        
        MBmodels_allVar{end+1} = MBmodel_varN_varW_varTh;
        allTestedOdors = PN;
				
    end

        % 
		% %% calculate KC responses (Y) to all odors
        % 
		% Y = calculateKCresponse(MBmodel,PNtrials);
        % totalKCactivity = 


