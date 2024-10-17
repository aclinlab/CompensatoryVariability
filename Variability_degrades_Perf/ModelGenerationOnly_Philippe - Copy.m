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
ll=5; % trials, each corresponding to different mushroom body instance
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
Sigmoid_derivative = @(x) exp(-0.1.*x)./((1.+exp(-0.1.*x)).^2)

C_SoftMax=1e-4;
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

PNs=[PNs_1 PNs_2_]; % now this is the 24 by 200 matrix (24 ORNs 200 odors)

% loops over different number of input odors, N=50,100,150,200.
% testing models at different number of odors, Fig2C1

% for mods=1:modss
mods = 2   
    odors= mulOd*(mods/modss);
    
    % input fictitious odors sampled from Hallem-Carlson
    x=PNs(:,1:odors);
    
    % original inputs from Hallem-Carlson
    %    x=PN(:,1:odors);
    
    % random fly instantiations
    for randomTrials=1:ll
        
        % random split of 'good' and 'bad' odors assignments: class 1 Vs.
        % class2
        classAction1=randsample([1:odors],round(odors/2),'false');
        
        APLgainP= zeros(1,6);
        APLgains=zeros(1,9);
        
        n =2000; % number of neurons in the hidden layer/ KCs
        
        m=24;  %number of dimensions in the input data/ PNs
        
        %select number of claws randomly 
        % -> for models depicted with green in Fig2
        clawsNo=(normrnd(6,1.7,[1 n])); 
        clawsNo(clawsNo<2)=2;
        clawsNo(clawsNo>11)=11;
        
        HomogClaws= ones([1,n])*6;
        
        PNsperKC = round(clawsNo.*ones(1,n)); %scale factor, simply =1 here
        
        HomogPNsperKC= HomogClaws;
        
        % set up individual KC connectivity, sample either variable or
        % homogenously in terms of PN-KC connections
        for i=1:n %iterate through KCs
            %draws from 1:24 (correspond to PN ids)
            PnToKc{i} = randsample(m, PNsperKC(i), true); 
            HomogPnToKc{i}= randsample(m, HomogPNsperKC(i), true);
        end % random initilaization of the weights
        
        %initialize the weights matrix between KC and PN
        thisW = zeros(m, n); %just reserves a variable, value irrelevant
        thisW_HomogModel=zeros(m,n); %value irrelevant
        
        %% 4 more models
        
        % variable n (=clawNo) &/or theta, fixed w.
        thisW_varN= zeros(m,n); %for model with var clawCount, fixed w
        
        % variable w but same n &/or theta,
        thisW_varWonly= zeros(m,n); %weights for fixed var N/theta, var W
        
        %instantiate weights for models with variable clawCount N
        % 2 flavors: var N & var w ; var N &fixed w
        for k=1:n
            for j=1:length(PnToKc{k})
                
                whichPN = PnToKc{k}(j);
                
                % pick random weight from a log normal distribution that
                % roughtly fits the Turner distribution
                thisWeight = exp(-0.0507+0.3527*randn(1)); %randn=Norm(0,1)
                % variable n but same weight values
                thisweight_same_var_n=1;
                
				%model with variable N (#claws) and weights
                thisW(whichPN, k) = thisW(whichPN, k) + thisWeight;
				%model with variable N but homogeneous weights
                thisW_varN(whichPN,k)= thisW_varN(whichPN,k)+ thisweight_same_var_n;
                
            end
        end
        
        %instantiate weights for models with homogeneous clawCount N        
        for k=1:n
            for j=1:length(HomogPnToKc{k})
                
                whichPN_homog= HomogPnToKc{k}(j);
                
				thisWeightHomo=1; %% homogenous equal unity weights connecting KCs to PNs.
                thisWeight = exp(-0.0507+0.3527*randn(1));
                
				%model with homogenous N (#claws) and homog. weights
                thisW_HomogModel(whichPN_homog,k)= thisWeightHomo+ thisW_HomogModel(whichPN_homog,k);
                %model with homogenous N (#claws) but variable weights
                thisW_varWonly(whichPN_homog,k)= thisW_varWonly (whichPN_homog,k) + thisWeight;
                
            end
        end
        
        
        % tune Ctheta & APL gains for realizing the experimental
        % sparsity levels: 10% and 20% (without inhibition/ APL feedback blocked).
        
        thetaS_Ftheta=10;
        
        step=1.0;
        
        % initalizing empty arrays for the calculations below
		% Ftheta refers to fixed theta models
		% -> I am confident these could be moved into next for-loop (for each noise)
		% --> actually they appear to never be used, so flagged for removal
        ActivationsDummy=zeros(n,odors*numtrainingSamples);
        ActivationsEqualizeddummy=zeros(n,odors*numtrainingSamples); %used to tune Ctheta and alpha
        ActivationsHomogenousdummy=zeros(n,odors*numtrainingSamples);
        ActivationsDummy_Ftheta=zeros(n,odors*numtrainingSamples);
        ActivationsEqualizeddummy_Ftheta=zeros(n,odors*numtrainingSamples);
        ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numtrainingSamples);
        Ydummy=zeros(n,odors*numtrainingSamples);
        YEqualizeddummy=zeros(n,odors*numtrainingSamples); %used to tune Ctheta and alpha
        YHomogdummy=zeros(n,odors*numtrainingSamples);
        
        Ydummy_Ftheta=zeros(n,odors*numtrainingSamples);
        YEqualizeddummy_Ftheta=zeros(n,odors*numtrainingSamples);
        YHomogdummy_Ftheta=zeros(n,odors*numtrainingSamples);
        
        %% flagging if the sparsity constraints are met: TRUE
        constraints=0;
        
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
                tic
                
                %start the parallel pool
                % spmd
                %
                %     if labindex==1
				%% begin optimisation procedure to get C_theta and APLgain values that satisfy coding level constraints
				%% for the variable weight, variable theta model.
                CLevelP=[];
                INHAbs_CLP=[];
                %thetaS is spike threshold in KCs (n=2000 KCs)
                T=10;
                thetaS=abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
                thetaS(thetaS>70)=70;
                thetaS(thetaS<0.01)=0.01;
				
				%n=2000 KCs, numtrainingSamples=15
                A=zeros(n,odors*numtrainingSamples); %A is excitation to KC_i from PNs to odor k
                Y=zeros(n,odors*numtrainingSamples); %resulting activity of KC without APL feedback
                eta=1; %scales adjustment steps for C_theta
                C_thetaS=1; %scaling factor to achieve the correct average coding level (10 resp. 20%)
				eta_2=0.0000001; %scales adjustment steps for ALPgain
                
                constraints=0; %%
                % while the Spar.constraints are NOT TRUE: repeat
                while(~constraints)
                    %calculate the KC responses to PN input for training trials
                    for trial = 1:(odors*numtrainingSamples)
                        
                        A(:,trial) = thisW'*PNtrials(:,trial); %excitation to KCs, first term of eq.3
						% calculate KC activity if at all (if A exceeds (adjusted) spiking threshold) (no inhibition)
                        Y(:,trial)=(( A(:,trial)-(C_thetaS.*thetaS) )>0 ).*( A(:,trial)-(C_thetaS.*thetaS));
                        codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n); %count proportion of active KCs
                    end
                    InhAbs_mSimp=mean(codingLevelDummy); %inhibition absent, 
                    
                    %% we want to change the theta so to achieve InhAbs_CL=2xCL
					% c.f. eq.14 and previous
					% CL stands for Coding Level, and we want 10% resp. 20% of KCs active with/without inhibition
                    depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2)); %shouldn't this be 0.8 instead?
                    depsi1_dy(isnan(depsi1_dy))=0;
                    depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(thetaS,1,odors*numtrainingSamples));
                    Grad= ((InhAbs_mSimp)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                    C_thetaS=C_thetaS - eta.*(Grad);
                    
                    if (C_thetaS<0)
                        error('the scale factor in the random model is -ve')
                    end
                    
                    %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                    % similar calculation, but including estimate of APLgain, and with adjusted C_theta
                    for trial = 1:(odors*numtrainingSamples)
                        ActivationsEqualizeddummy(:,trial) = thisW'*PNtrials(:,trial);
                        YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-...
							APLgainP(1)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-...
							(C_thetaS.*thetaS) )>0 ) .* ...
						  ( ActivationsEqualizeddummy(:,trial)- ...
						    APLgainP(1)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-...
							(C_thetaS.*thetaS) );
                        codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                    end
                    CLRand=mean(codingLevelEqualizedDummy);
                    
					%c.f. eq.18 and previous, adjusting APLgain alpha
                    dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                    dsig_dy(isnan(dsig_dy))=0;
                    dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                    dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                    Grad_alpha= ((CLRand)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                    APLgainP(1)= APLgainP(1)- eta_2*(Grad_alpha);
                    
                    %% checking if the constraints are met:
					% repeat calculations done during optimisation but with updated values C_theta, APLgain
                    %  first for zero-inhibition condition
					for trial = 1:(odors*numtrainingSamples)
                        Activations_S(:,trial) = thisW'*PNtrials(:,trial);
                        Y_S(:,trial)=(( Activations_S(:,trial)-(C_thetaS.*thetaS) )>0 )...
							.*( Activations_S(:,trial)-(C_thetaS.*thetaS) );
                        codingLevelDummy(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                    end
                    InhAbs_CL=mean(codingLevelDummy);
                    
					%  then for KC activity with APL inhibition
                    for trial = 1:(odors*numtrainingSamples)
                        Activation(:,trial) = thisW'*PNtrials(:,trial);
                        Y(:,trial)=(( Activation(:,trial)-(APLgainP(1))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS) )>0 ).*( Activation(:,trial)-APLgainP(1)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS) );
                        codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                    end
                    CL_=mean(codingLevelDummy);
                    
                    constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( abs(CL_-0.10)<0.01 );
                    
                end
				% end optimisation, store results
                CLevelP(end+1)=CL_; %appends to vector, empty until now (renewed for each noise level)
                thetaS=(C_thetaS.*thetaS) ;
                INHAbs_CLP(end+1)=InhAbs_CL;
             
			 toc
			 
			%%---------------------------------------------------------------------------%
				% for the homogenous weights, variable theta model.
				%for comments, refer to previous optimisation for variable weight model
                % Th=10;
                % thetaH=abs(normrnd(Th,Th*(5.6/21.5),[n 1])); %% avoid negative values of theta
                % thetaH(thetaH>70)=70;
                % thetaH(thetaH<0.01)=0.01;
                % constraints=0; %%
                % A=zeros(n,odors*numtrainingSamples);
                % Y=zeros(n,odors*numtrainingSamples);
                % eta=1;
                % C_thetaH=1;
                
                % while(~constraints)
                    
                    % for trial = 1:(odors*numtrainingSamples)
                        % A(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                        % Y(:,trial)=(( A(:,trial)- (C_thetaH*thetaH) )>0 ).*( A(:,trial)-(C_thetaH*thetaH));
                        % codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                    % end
                    % InhAbs_mHomog=mean(codingLevelDummy);
                    % we want to change the theta so to achieve InhAbs_CL=2xCL
                    % depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                    % depsi1_dy(isnan(depsi1_dy))=0;
                    % depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaH,1,odors*numtrainingSamples));
                    % Grad= ((InhAbs_mHomog)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                    % C_thetaH=C_thetaH - eta.*(Grad);
                    
                    % if (C_thetaH<0)
                        % error('the scale factor in black model is -ve')
                    % end
                    
                    % allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                    % eta_2=0.0000001;
                    % for trial = 1:(odors*numtrainingSamples)
                        % ActivationsEqualizeddummy(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                        % YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(2))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)- (C_thetaH*thetaH) )>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(2)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaH*thetaH));
                        % codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                    % end
                    % CLHomogVarTh=mean(codingLevelEqualizedDummy);
                    % dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                    % dsig_dy(isnan(dsig_dy))=0;
                    % dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                    % dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                    % Grad_alpha= ((CLHomogVarTh)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                    % APLgainP(2)= APLgainP(2)- eta_2*(Grad_alpha);
                    
                    % check if the constraints are satisfied
                    % for trial = 1:(odors*numtrainingSamples)
                        
                        % Activations_H(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                        % Y_H(:,trial)=(( Activations_H(:,trial)-(C_thetaH*thetaH))>0 ).*( Activations_H(:,trial)-(C_thetaH*thetaH));
                        % codingLevelDummy(trial)=  (sum(Y_H(:,trial)>0,1)/n);
                    % end
                    % InhAbs_CL=mean(codingLevelDummy);
                    
                    % for trial = 1:(odors*numtrainingSamples)
                        % Activation(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                        % Y(:,trial)=(( Activation(:,trial)-(APLgainP(2))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH*thetaH))>0 ).*( Activation(:,trial)-APLgainP(2)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH*thetaH));
                        % codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                    % end
                    % CL_=mean(codingLevelDummy);
                    
                    % constraints= (abs( (InhAbs_CL/CL_) - 2.0)<0.2) &( abs(CL_-0.10)<0.01 );
                % end
                %appends results to same vector
                % CLevelP(end+1)=CL_;
                % INHAbs_CLP(end+1)=InhAbs_CL;
                % thetaH=(C_thetaH*thetaH);
                
                
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
                
				
				%% calculate KC responses (Y) to all odors
				
				% Activation variable refers to the first term in eq.3, i.e. the excitatory drive to KCs before inhibition or thresholding
                % This is done first because the spiking activity combines the exc. and the APL inhib., whihc depends on total exc. as well
				Activations=zeros(n,odors*numTrials);
                % ActivationsHomog=zeros(n,odors*numTrials);
                % Activations_varW_fixedN= zeros(n,odors*numTrials);
                % Activations_varN_fixedW= zeros(n,odors*numTrials);
                % Activations_Ftheta=zeros(n,odors*numTrials);
                % ActivationsHomog_Ftheta=zeros(n,odors*numTrials);
                % Activations_varW_fixedN_and_theta= zeros(n,odors*numTrials);
                % Activations_varN_fixedW_and_theta= zeros(n,odors*numTrials);
                
                % Y refers to the spiking/output activity of KCs, incorporating APL inhibition and spiking thresholds
                Y=zeros(n,odors*numTrials); %model 1
                % YHomog=zeros(n,odors*numTrials); %model 2
                % Y_varw_fixedN= zeros(n,odors*numTrials); %model 3
                % Y_varN_fixedW= zeros(n,odors*numTrials); 
                % Y_Ftheta=zeros(n,odors*numTrials); 
                % YHomog_Ftheta=zeros(n,odors*numTrials); %model 4
                % Y_varw_fixedN_and_theta= zeros(n,odors*numTrials);
                % Y_varN_fixedW_and_theta= zeros(n,odors*numTrials);
                
				% double-checked all variable names with parameter combinations -> OK
                for trial = 1:(odors*numTrials)
                    
                    %ActivationsHomog_Ftheta(:,trial) = thisW_HomogModel'*PNtrials(:,trial  );
                    %YHomog_Ftheta(:,trial)=(( ActivationsHomog_Ftheta(:,trial)-(APLgains(4))*repmat(sum(ActivationsHomog_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomog_Ftheta(:,trial)-APLgains(4)*repmat(sum(ActivationsHomog_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);
                    
                    % Activations_Ftheta(:,trial) = thisW'*PNtrials(:,trial );
                    % Y_Ftheta(:,trial)=(( Activations_Ftheta(:,trial)-(APLgains(3))*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta)>0 ).*( Activations_Ftheta(:,trial)-APLgains(3)*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta);
                    
                    % Activations_varW_fixedN_and_theta(:,trial) = thisW_varWonly'*PNtrials(:,trial );
                    % Y_varw_fixedN_and_theta(:,trial)=(( Activations_varW_fixedN_and_theta(:,trial)-(APLgains(6) )*repmat(sum(Activations_varW_fixedN_and_theta(:,trial),1),n,1)-theta_varw_fixedN_scalar)>0 ).*( Activations_varW_fixedN_and_theta(:,trial)-APLgains(6)*repmat(sum(Activations_varW_fixedN_and_theta(:,trial),1),n,1)-theta_varw_fixedN_scalar);
                    
                    % Activations_varN_fixedW_and_theta(:,trial) = thisW_varN'*PNtrials(:,trial );
                    % Y_varN_fixedW_and_theta(:,trial)=(( Activations_varN_fixedW_and_theta(:,trial)-(APLgains(8) )*repmat(sum(Activations_varN_fixedW_and_theta(:,trial),1),n,1)-theta_varN_fixedW_scalar)>0 ).*( Activations_varN_fixedW_and_theta(:,trial)-APLgains(8)*repmat(sum(Activations_varN_fixedW_and_theta(:,trial),1),n,1)-theta_varN_fixedW_scalar);
                    
                    % ActivationsHomog(:,trial) = thisW_HomogModel'*PNtrials(:,trial );
                    % YHomog(:,trial)=(( ActivationsHomog(:,trial)-(APLgains(2))*repmat(sum(ActivationsHomog(:,trial),1),n,1)-thetaH)>0 ).*( ActivationsHomog(:,trial)-APLgains(2)*repmat(sum(ActivationsHomog(:,trial),1),n,1)-thetaH);
                    
                    Activations(:,trial) = thisW'*PNtrials(:,trial );
                    Y(:,trial)=(( Activations(:,trial)-(APLgains(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgains(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
                    
                    % Activations_varW_fixedN(:,trial) = thisW_varWonly'*PNtrials(:,trial );
                    % Y_varw_fixedN(:,trial)=(( Activations_varW_fixedN(:,trial)-(APLgains(5) )*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN)>0 ).*( Activations_varW_fixedN(:,trial)-APLgains(5)*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN);
                    
                    % Activations_varN_fixedW(:,trial) = thisW_varN'*PNtrials(:,trial );
                    % Y_varN_fixedW(:,trial)=(( Activations_varN_fixedW(:,trial)-(APLgains(7) )*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw)>0 ).*( Activations_varN_fixedW(:,trial)-APLgains(7)*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw);
                    
                end
    end
    
%------------------------------------------------END OF CODE



