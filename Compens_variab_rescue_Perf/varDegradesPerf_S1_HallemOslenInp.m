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
ll=30;
numtrainingSamples=15;

% at 4 scales of 200 odors, in increment of 50: 50,100,150,200.
% to test at just one scale (like 100 odors): simply change
% modss=1, and mulOd=100
modss=1;
mulOd=110;
odors=mulOd;

numTrials = 30; %number of noisy trials per each odors
PN = hallem_olsen(1:110,:)';
k=0;

C_SoftMax=1;
% levels of noise added to inputs
NScales=1;

% learning rates for the KC-MBON weights. 7 scales
lrs=7;

% scales for indetermincy in the decision making (Softmax function--> Step
% function) 5 scales of C investigated in Fig2C4
Crange=1;


% loops over different number of input odors, N=50,100,150,200.
% testing models at different number of odors, Fig2C1
for mods=1:modss
    
    odors= mulOd*(mods/modss);
    
    %original inputs from Hallem-Carlson
    x=PN(:,1:odors);
    
    % random fly instantiations
    for randomTrials=1:ll
        randomTrials
        % load the 100 odors from 2E
        
        
        % random split of 'good' and 'bad' odors assignments: class 1 Vs.
        % class2
        classAction1=randsample([1:odors],round(odors/2),'false');
        
        APLgainP= zeros(1,6);
        APLgains=zeros(1,9);
        
        n =2000; % number of neurons in the hidden layer/ KCs
        
        m=24;  %number of dimensions in the input data/ PNs
        
        clawsNo=(normrnd(6,1.7,[1 n])); %select number of claws randomly
        clawsNo(clawsNo<2)=2;
        clawsNo(clawsNo>11)=11;
        
        HomogClaws= ones([1,n])*6;
        
        PNsperKC = round(clawsNo.*ones(1,n));
        
        HomogPNsperKC= HomogClaws;
        
        for i=1:n
            PnToKc{i} = randsample(m, PNsperKC(i), true);
            HomogPnToKc{i}= randsample(m, HomogPNsperKC(i), true);
        end % random initilaization of the weights
        
        %initialize the weights matrix between KC and PN
        thisW = zeros(m, n);
        thisW_HomogModel=zeros(m,n);
        
        %% 4 more models
        
        % variable n &/or theta, fixed w.
        thisW_varN= zeros(m,n);
        
        % variable w but same n &/or theta,
        thisW_varWonly= zeros(m,n);
        
        
        for k=1:n
            
            for j=1:length(PnToKc{k})
                
                whichPN = PnToKc{k}(j);
                
                % pick random weight from a log normal distribution that
                % roughtly fits the Turner distribution
                
                thisWeight = exp(-0.0507+0.3527*randn(1));
                
                % variable n but same weight values
                thisweight_same_var_n=1;
                
                
                thisW(whichPN, k) = thisW(whichPN, k) + thisWeight;
                thisW_varN(whichPN,k)= thisW_varN(whichPN,k)+ thisweight_same_var_n;
                
            end
        end
        
        
        for k=1:n
            
            for j=1:length(HomogPnToKc{k})
                
                whichPN_homog= HomogPnToKc{k}(j);
                
                thisWeightHomo=1; %% homogenous equal unity weights connecting KCs to PNs.
                
                thisWeight = exp(-0.0507+0.3527*randn(1));
                
                thisW_HomogModel(whichPN_homog,k)= thisWeightHomo+ thisW_HomogModel(whichPN_homog,k);
                
                thisW_varWonly(whichPN_homog,k)= thisW_varWonly (whichPN_homog,k) + thisWeight;
                
            end
        end
        
        
        % tune Ctheta & APL gains for realizing the experimental
        % sparsity levels: 10% and 20% (without inhibition/ APL feedback blocked).
        
        thetaS_Ftheta=10;
        
        step=1.0;
        
        % initalizing empty arrays for the calculations below
        ActivationsDummy=zeros(n,odors*numtrainingSamples);
        ActivationsEqualizeddummy=zeros(n,odors*numtrainingSamples);
        ActivationsHomogenousdummy=zeros(n,odors*numtrainingSamples);
        ActivationsDummy_Ftheta=zeros(n,odors*numtrainingSamples);
        ActivationsEqualizeddummy_Ftheta=zeros(n,odors*numtrainingSamples);
        ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numtrainingSamples);
        Ydummy=zeros(n,odors*numtrainingSamples);
        YEqualizeddummy=zeros(n,odors*numtrainingSamples);
        YHomogdummy=zeros(n,odors*numtrainingSamples);
        
        Ydummy_Ftheta=zeros(n,odors*numtrainingSamples);
        YEqualizeddummy_Ftheta=zeros(n,odors*numtrainingSamples);
        YHomogdummy_Ftheta=zeros(n,odors*numtrainingSamples);
        
        %% flagging if the sparsity constraints are met: TRUE
        constraints=0;
        
        % testing models at different scales of input noise: Fig2C3
        for noiseScale=1:NScales
            
            noise=1; %((noiseScale^2)/10);
            PNtrials = zeros(24, odors, numTrials);
            SNR=zeros(24,odors,numTrials);
            PNtrials(:,:,1) =x;
            for t = 1:numTrials-1
                PNtrials(:,:,t+1) = x + ...
                    getPNStdevBhandawat(x) .* ...
                    noise.*randn(24, odors);
                SNR(:,:,t)= (x./(getPNStdevBhandawat(x) .* ...
                    noise.*randn(24, odors))).^2;
                SNR(:,:,t)=log(SNR(:,:,t));
            end
            
            SNR_averageAcrossTrials_andPNs(mods,randomTrials,noiseScale)= mean(SNR(:));
            
            % making sure that PN responses are +ve
            PNtrials(PNtrials<0)=0;
            
            %rescaling the PN responses in range from 0 to 5
            PNtrials=rescale(PNtrials,0,5);
            tic
            
            %start the parallel pool
            % spmd
            %
            %     if labindex==1
            CLevelP=[];
            INHAbs_CLP=[];
            T=10;
            thetaS=abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            thetaS(thetaS>70)=70;
            thetaS(thetaS<0.01)=0.01;
            constraints=0; %%
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            eta=1;
            C_thetaS=1;
            
            % while the Spar.constraints are NOT TRUE: repeat
            while(~constraints)
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    A(:,trial) = thisW'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)-(C_thetaS.*thetaS) )>0 ).*( A(:,trial)-(C_thetaS.*thetaS));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                InhAbs_mSimp=mean(codingLevelDummy);
                
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(thetaS,1,odors*numtrainingSamples));
                Grad= ((InhAbs_mSimp)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                C_thetaS=C_thetaS - eta.*(Grad);
                
                if (C_thetaS<0)
                    error('the scale factor in the random model is -ve')
                end
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsEqualizeddummy(:,trial) = thisW'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(1))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)- (C_thetaS.*thetaS) )>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(1)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaS.*thetaS) );
                    codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                end
                CLRand=mean(codingLevelEqualizedDummy);
                
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CLRand)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgainP(1)= APLgainP(1)- eta_2*(Grad_alpha);
                
                %% checking if the constraints are met:
                for trial = 1:(odors*numtrainingSamples)
                    Activations_S(:,trial) = thisW'*PNtrials(:,trial);
                    Y_S(:,trial)=(( Activations_S(:,trial)-(C_thetaS.*thetaS) )>0 ).*( Activations_S(:,trial)-(C_thetaS.*thetaS) );
                    codingLevelDummy(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgainP(1))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS) )>0 ).*( Activation(:,trial)-APLgainP(1)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS) );
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevelDummy);
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( abs(CL_-0.10)<0.01 );
                
                
            end
            CLevelP(end+1)=CL_;
            thetaS=(C_thetaS.*thetaS) ;
            INHAbs_CLP(end+1)=InhAbs_CL;
            
            %        elseif labindex==2
            % for the homogenous weights, variable theta model.
            Th=10;
            thetaH=abs(normrnd(Th,Th*(5.6/21.5),[n 1])); %% avoid negative values of theta
            thetaH(thetaH>70)=70;
            thetaH(thetaH<0.01)=0.01;
            constraints=0; %%
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            eta=1;
            C_thetaH=1;
            
            while(~constraints)
                
                for trial = 1:(odors*numtrainingSamples)
                    A(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)- (C_thetaH*thetaH) )>0 ).*( A(:,trial)-(C_thetaH*thetaH));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                InhAbs_mHomog=mean(codingLevelDummy);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaH,1,odors*numtrainingSamples));
                Grad= ((InhAbs_mHomog)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                C_thetaH=C_thetaH - eta.*(Grad);
                
                if (C_thetaH<0)
                    error('the scale factor in black model is -ve')
                end
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsEqualizeddummy(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(2))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)- (C_thetaH*thetaH) )>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(2)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaH*thetaH));
                    codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                end
                CLHomogVarTh=mean(codingLevelEqualizedDummy);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CLHomogVarTh)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgainP(2)= APLgainP(2)- eta_2*(Grad_alpha);
                
                %% check if the constraints are satisfied
                for trial = 1:(odors*numtrainingSamples)
                    
                    Activations_H(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y_H(:,trial)=(( Activations_H(:,trial)-(C_thetaH*thetaH))>0 ).*( Activations_H(:,trial)-(C_thetaH*thetaH));
                    codingLevelDummy(trial)=  (sum(Y_H(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgainP(2))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH*thetaH))>0 ).*( Activation(:,trial)-APLgainP(2)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH*thetaH));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevelDummy);
                
                constraints= (abs( (InhAbs_CL/CL_) - 2.0)<0.2) &( abs(CL_-0.10)<0.01 );
            end
            
            CLevelP(end+1)=CL_;
            INHAbs_CLP(end+1)=InhAbs_CL;
            thetaH=(C_thetaH*thetaH);
            
            %% ---------------------------------------------------------------------
            %%%%%%%  Optimization  FOR THE FIXED THREHSOLD MODELS %%%%%%
            
            %     elseif labindex==3
            % Random Weights fixed thresholds model
            constraints=0;
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            C_thetaS_Ftheta=1;
            while(~constraints)
                
                eta=1;
                for trial = 1:(odors*numtrainingSamples)
                    A(:,trial) = thisW'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)- (C_thetaS_Ftheta*thetaS_Ftheta) )>0 ).*( A(:,trial)-(C_thetaS_Ftheta*thetaS_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                InhAbs_mSimp=mean(codingLevelDummy);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaS_Ftheta,n,odors*numtrainingSamples));
                Grad= ((InhAbs_mSimp)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                C_thetaS_Ftheta=C_thetaS_Ftheta - eta.*(Grad);
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsEqualizeddummy(:,trial) = thisW'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(3))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(3)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta));
                    codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                end
                CLRandFixedTh=mean(codingLevelEqualizedDummy);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CLRandFixedTh)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgainP(3)= APLgainP(3)- eta_2*(Grad_alpha);
                
                %% check if the constraints are satisfied
                for trial = 1:(odors*numtrainingSamples)
                    Activations_SF(:,trial) = thisW'*PNtrials(:,trial);
                    Y_SF(:,trial)=(( Activations_SF(:,trial)-(C_thetaS_Ftheta*thetaS_Ftheta))>0 ).*( Activations_SF(:,trial)-(C_thetaS_Ftheta*thetaS_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y_SF(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgainP(3))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta))>0 ).*( Activation(:,trial)-APLgainP(3)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevelDummy);
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01 );
                
            end
            
            CLevelP(end+1)=CL_;
            INHAbs_CLP(end+1)=InhAbs_CL;
            thetaS_Ftheta=(C_thetaS_Ftheta*thetaS_Ftheta);
            
            
            %        elseif labindex==4
            % Fixed weights &thresh. model
            constraints=0;
            thetaH_Ftheta=10+rand(1);
            C_thetaHF=1;
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            HomoFtheta_Inh=0;
            
            while(~constraints)
                
                eta=1;
                for trial = 1:(odors*numtrainingSamples)
                    A(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)- (C_thetaHF*thetaH_Ftheta) )>0 ).*( A(:,trial)-(C_thetaHF*thetaH_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                InhAbs_mHomog=mean(codingLevelDummy);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*repmat(thetaH_Ftheta,n,odors*numtrainingSamples);
                Grad= ((InhAbs_mHomog)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                C_thetaHF=C_thetaHF - eta.*(Grad);
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsEqualizeddummy(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(4))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaHF*thetaH_Ftheta))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(4)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaHF*thetaH_Ftheta));
                    codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                end
                CLAllFixed=mean(codingLevelEqualizedDummy);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CLAllFixed)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgainP(4)= APLgainP(4)- eta_2*(Grad_alpha);
                
                %% check if the constraints are satisfied
                for trial = 1:(odors*numtrainingSamples)
                    Activations_H(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y_H(:,trial)=(( Activations_H(:,trial)-(C_thetaHF*thetaH_Ftheta))>0 ).*( Activations_H(:,trial)-(C_thetaHF*thetaH_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y_H(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgainP(4))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaHF*thetaH_Ftheta))>0 ).*( Activation(:,trial)-APLgainP(4)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaHF*thetaH_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevelDummy);
                
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2) &( abs(CL_-0.10)<0.01 );
                
                
                if(abs(InhAbs_CL-0.37)<0.000001)
                    break;
                end
                
            end
            CLevelP(end+1)=CL_;
            INHAbs_CLP(end+1)=InhAbs_CL;
            thetaH_Ftheta=(C_thetaHF*thetaH_Ftheta);
            
            %     end [CLevelP{1},CLevelP{2},CLevelP{3},CLevelP{4}];
            % end [INHAbs_CLP{1},INHAbs_CLP{2},INHAbs_CLP{3},INHAbs_CLP{4}];
            Clevels=CLevelP;
            INHAbs_CL=INHAbs_CLP;
            
            for i=1:4
                APLgains(i)=APLgainP(i);
                %                     = temp(i);
            end
            
            %                 thetaH=thetaH{2};
            %                 thetaS=thetaS{1};
            %                 thetaS_Ftheta=thetaS_Ftheta{3};
            %                 thetaH_Ftheta=thetaH_Ftheta{4};
            %
            
            %% ---  % Model 5, same n, variable W & theta
            
            constraints=0; %% again for the variable theta, random and homog models
            theta_varw_fixedN= abs(normrnd(12,12*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_varw_fixedN(theta_varw_fixedN>70)=70;
            theta_varw_fixedN(theta_varw_fixedN<0.01)=0.01;
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            Activation=zeros(n,odors*numtrainingSamples);
            Activations_S=zeros(n,odors*numtrainingSamples);
            Y_S=zeros(n,odors*numtrainingSamples);
            ActivationsEqualizeddummy=zeros(n,odors*numtrainingSamples);
            YEqualizeddummy=zeros(n,odors*numtrainingSamples);
            C_theta5=1;
            eta=5;
            
            while(~constraints)
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    A(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)- (C_theta5.*theta_varw_fixedN) )>0 ).*( A(:,trial)-(C_theta5.*theta_varw_fixedN));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                    
                end
                InhAbs_fixedN_varw=mean(codingLevel);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                
                depsi1_dtheta= -(Y>0).* depsi1_dy.*repmat(theta_varw_fixedN,1,odors*numtrainingSamples);
                
                
                Grad= ((InhAbs_fixedN_varw)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                
                C_theta5=C_theta5 - eta.*(Grad);
                
                if(C_theta5<0)
                    error('Ctheta in model 5 cannot be negative')
                end
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                
                eta_2=0.0000001;
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    ActivationsEqualizeddummy(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgains(5))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta5.*theta_varw_fixedN))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgains(5)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta5.*theta_varw_fixedN));
                    codingLevelEqualized(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                    
                end
                CL_varw_fixedN=mean(codingLevelEqualized);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_varw_fixedN)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(5)= APLgains(5)- eta_2*(Grad_alpha);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activations_S(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    Y_S(:,trial)=(( Activations_S(:,trial)-(C_theta5.*theta_varw_fixedN))>0 ).*( Activations_S(:,trial)-(C_theta5.*theta_varw_fixedN));
                    codingLevel(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevel);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgains(5))*repmat(sum(Activation(:,trial),1),n,1)-(C_theta5.*theta_varw_fixedN))>0 ).*( Activation(:,trial)-APLgains(5)*repmat(sum(Activation(:,trial),1),n,1)-(C_theta5.*theta_varw_fixedN));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevel);
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( abs(CL_-0.10)<0.01 );
                
            end
            
            Clevels(5)=CL_;
            INHAbs_CL(5)=InhAbs_CL;
            theta_varw_fixedN=(C_theta5.*theta_varw_fixedN);
            
            %% ------- % Model 6, same n &theta, variable W ONLY
            constraints=0;
            theta_varw_fixedN_scalar= 12+ rand(1);
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            C_theta6=1;
            
            eta=0.1;
            
            while(~constraints)
                
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    A(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)- (C_theta6*theta_varw_fixedN_scalar) )>0 ).*( A(:,trial)-(C_theta6*theta_varw_fixedN_scalar));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                    
                end
                InhAbs_fixedN_varwS=mean(codingLevel);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                
                depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(theta_varw_fixedN_scalar,n,odors*numtrainingSamples));
                Grad= ((InhAbs_fixedN_varwS)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                C_theta6=C_theta6 - eta.*(Grad);
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                for trial = 1:(odors*numtrainingSamples)
                    
                    ActivationsEqualizeddummy(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgains(6))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta6*theta_varw_fixedN_scalar))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgains(6)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta6*theta_varw_fixedN_scalar));
                    codingLevelEqualized(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                    
                end
                CL_varw_fixedNS=mean(codingLevelEqualized);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_varw_fixedNS)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(6)= APLgains(6)- eta_2*(Grad_alpha);
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    Activations_S(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    Y_S(:,trial)=(( Activations_S(:,trial)-(C_theta6*theta_varw_fixedN_scalar))>0 ).*( Activations_S(:,trial)-(C_theta6*theta_varw_fixedN_scalar));
                    codingLevel(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevel);
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    Activation(:,trial) = thisW_varWonly'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgains(6))*repmat(sum(Activation(:,trial),1),n,1)-(C_theta6*theta_varw_fixedN_scalar))>0 ).*( Activation(:,trial)-APLgains(6)*repmat(sum(Activation(:,trial),1),n,1)-(C_theta6*theta_varw_fixedN_scalar));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevel);
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( abs(CL_-0.10)<0.01 );
                
            end
            
            Clevels(6)=CL_;
            INHAbs_CL(6)=InhAbs_CL;
            theta_varw_fixedN_scalar=(C_theta6*theta_varw_fixedN_scalar);
            
            %% ------ % Model 7, same W, variable N & theta
            
            constraints=0;
            theta_varN_fixedw= abs(normrnd(12,12*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_varN_fixedw(theta_varN_fixedw>70)=70;
            theta_varN_fixedw(theta_varN_fixedw<0.01)=0.01;
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            C_theta7=1;
            eta=1;
            while(~constraints)
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    A(:,trial) = thisW_varN'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)- (C_theta7*theta_varN_fixedw) )>0 ).*( A(:,trial)-(C_theta7*theta_varN_fixedw));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                    
                end
                InhAbss=mean(codingLevel);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(theta_varN_fixedw,1,odors*numtrainingSamples));
                Grad= ((InhAbss)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                C_theta7=C_theta7 - eta.*(Grad);
                if (C_theta7<0)
                    error('c theta in model 7 cannot be negative')
                end
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsEqualizeddummy(:,trial) = thisW_varN'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgains(7))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta7*theta_varN_fixedw))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgains(7)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta7*theta_varN_fixedw));
                    codingLevelEqualized(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                end
                CL_varN=mean(codingLevelEqualized);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_varN)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(7)= APLgains(7)- eta_2*(Grad_alpha);
                
                
                for trial = 1:(odors*numtrainingSamples)
                    Activations_S(:,trial) = thisW_varN'*PNtrials(:,trial);
                    Y_S(:,trial)=(( Activations_S(:,trial)-(C_theta7*theta_varN_fixedw))>0 ).*( Activations_S(:,trial)-(C_theta7*theta_varN_fixedw));
                    codingLevel(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevel);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW_varN'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgains(7))*repmat(sum(Activation(:,trial),1),n,1)-(C_theta7*theta_varN_fixedw))>0 ).*( Activation(:,trial)-APLgains(7)*repmat(sum(Activation(:,trial),1),n,1)-(C_theta7*theta_varN_fixedw));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevel);
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( abs(CL_-0.10)<0.01 );
                
            end
            
            Clevels(7)=CL_;
            INHAbs_CL(7)=InhAbs_CL;
            theta_varN_fixedw=(C_theta7*theta_varN_fixedw);
            
            %% ----- % Model 8, same W and theta, variable N
            
            constraints=0;
            theta_varN_fixedW_scalar= 11+ rand(1);
            A=zeros(n,odors*numtrainingSamples);
            Y=zeros(n,odors*numtrainingSamples);
            C_theta8=1;
            eta=1;
            
            while(~constraints)
                
                for trial = 1:(odors*numtrainingSamples)
                    A(:,trial) = thisW_varN'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)-  (C_theta8*theta_varN_fixedW_scalar) )>0 ).*( A(:,trial)- (C_theta8*theta_varN_fixedW_scalar));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                InhAbss=mean(codingLevel);
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(theta_varN_fixedW_scalar,n,odors*numtrainingSamples));
                Grad= ((InhAbss)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                C_theta8=C_theta8 - eta.*(Grad);
                
                if(C_theta8<0)
                    error('Ctheta 8 cant be -ve')
                end
                
                %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
                eta_2=0.0000001;
                
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsEqualizeddummy(:,trial) = thisW_varN'*PNtrials(:,trial);
                    YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgains(8))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)- (C_theta8*theta_varN_fixedW_scalar))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgains(8)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)- (C_theta8*theta_varN_fixedW_scalar));
                    codingLevelEqualized(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                end
                CL_varNS=mean(codingLevelEqualized);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_varNS)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(8)= APLgains(8)- eta_2*(Grad_alpha);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activations_S(:,trial) = thisW_varN'*PNtrials(:,trial);
                    Y_S(:,trial)=(( Activations_S(:,trial)- (C_theta8*theta_varN_fixedW_scalar))>0 ).*( Activations_S(:,trial)- (C_theta8*theta_varN_fixedW_scalar));
                    codingLevel(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevel);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activation(:,trial) = thisW_varN'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgains(8))*repmat(sum(Activation(:,trial),1),n,1)- (C_theta8*theta_varN_fixedW_scalar))>0 ).*( Activation(:,trial)-APLgains(8)*repmat(sum(Activation(:,trial),1),n,1)- (C_theta8*theta_varN_fixedW_scalar));
                    codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevel);
                
                constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( abs(CL_-0.10)<0.01 );
                
            end
            Clevels(8)=CL_;
            INHAbs_CL(8)=InhAbs_CL;
            theta_varN_fixedW_scalar= (C_theta8*theta_varN_fixedW_scalar);
            
            toc
            
            %% if the loops broke without the constraints being statisfied. then
            %% print an error message
            
            if(any(APLgains<0) || any([theta_varN_fixedW_scalar,theta_varw_fixedN_scalar,thetaH_Ftheta,thetaS_Ftheta]<0) )
                
                error( 'constraints cant be statisfied ' )
            end
            
            %% save the parameters for each fly, connectivity weights, spiking thresholds, inhibition Gains,for all 8 models
            
            if(odors==110)
                save(strcat(['VarDegradesPerf_fly_wNoise_HOInp',num2str(randomTrials),num2str(noiseScale)]) , 'thisW_HomogModel','APLgains',...
                    'thetaH_Ftheta', 'thisW','thetaS_Ftheta' ,'thisW_varWonly','thisW_varN','theta_varw_fixedN_scalar','theta_varN_fixedW_scalar',...
                    'thetaH', 'thetaS','theta_varw_fixedN', 'theta_varN_fixedw','PNtrials','classAction1' );
                
            end
            
            Activations=zeros(n,odors*numTrials);
            ActivationsHomog=zeros(n,odors*numTrials);
            Activations_varW_fixedN= zeros(n,odors*numTrials);
            Activations_varN_fixedW= zeros(n,odors*numTrials);
            Activations_Ftheta=zeros(n,odors*numTrials);
            ActivationsHomog_Ftheta=zeros(n,odors*numTrials);
            Activations_varW_fixedN_and_theta= zeros(n,odors*numTrials);
            Activations_varN_fixedW_and_theta= zeros(n,odors*numTrials);
            
            
            Y=zeros(n,odors*numTrials);
            YHomog=zeros(n,odors*numTrials);
            Y_varw_fixedN= zeros(n,odors*numTrials);
            Y_varN_fixedW= zeros(n,odors*numTrials);
            Y_Ftheta=zeros(n,odors*numTrials);
            YHomog_Ftheta=zeros(n,odors*numTrials);
            Y_varw_fixedN_and_theta= zeros(n,odors*numTrials);
            Y_varN_fixedW_and_theta= zeros(n,odors*numTrials);
            
            
            for trial = 1:(odors*numTrials)
                
                ActivationsHomog_Ftheta(:,trial) = thisW_HomogModel'*PNtrials(:,trial  );
                YHomog_Ftheta(:,trial)=(( ActivationsHomog_Ftheta(:,trial)-(APLgains(4))*repmat(sum(ActivationsHomog_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomog_Ftheta(:,trial)-APLgains(4)*repmat(sum(ActivationsHomog_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);
                
                Activations_Ftheta(:,trial) = thisW'*PNtrials(:,trial );
                Y_Ftheta(:,trial)=(( Activations_Ftheta(:,trial)-(APLgains(3))*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta)>0 ).*( Activations_Ftheta(:,trial)-APLgains(3)*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta);
                
                Activations_varW_fixedN_and_theta(:,trial) = thisW_varWonly'*PNtrials(:,trial );
                Y_varw_fixedN_and_theta(:,trial)=(( Activations_varW_fixedN_and_theta(:,trial)-(APLgains(6) )*repmat(sum(Activations_varW_fixedN_and_theta(:,trial),1),n,1)-theta_varw_fixedN_scalar)>0 ).*( Activations_varW_fixedN_and_theta(:,trial)-APLgains(6)*repmat(sum(Activations_varW_fixedN_and_theta(:,trial),1),n,1)-theta_varw_fixedN_scalar);
                
                Activations_varN_fixedW_and_theta(:,trial) = thisW_varN'*PNtrials(:,trial );
                Y_varN_fixedW_and_theta(:,trial)=(( Activations_varN_fixedW_and_theta(:,trial)-(APLgains(8) )*repmat(sum(Activations_varN_fixedW_and_theta(:,trial),1),n,1)-theta_varN_fixedW_scalar)>0 ).*( Activations_varN_fixedW_and_theta(:,trial)-APLgains(8)*repmat(sum(Activations_varN_fixedW_and_theta(:,trial),1),n,1)-theta_varN_fixedW_scalar);
                
                ActivationsHomog(:,trial) = thisW_HomogModel'*PNtrials(:,trial );
                YHomog(:,trial)=(( ActivationsHomog(:,trial)-(APLgains(2))*repmat(sum(ActivationsHomog(:,trial),1),n,1)-thetaH)>0 ).*( ActivationsHomog(:,trial)-APLgains(2)*repmat(sum(ActivationsHomog(:,trial),1),n,1)-thetaH);
                
                Activations(:,trial) = thisW'*PNtrials(:,trial );
                Y(:,trial)=(( Activations(:,trial)-(APLgains(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgains(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
                
                Activations_varW_fixedN(:,trial) = thisW_varWonly'*PNtrials(:,trial );
                Y_varw_fixedN(:,trial)=(( Activations_varW_fixedN(:,trial)-(APLgains(5) )*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN)>0 ).*( Activations_varW_fixedN(:,trial)-APLgains(5)*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN);
                
                Activations_varN_fixedW(:,trial) = thisW_varN'*PNtrials(:,trial );
                Y_varN_fixedW(:,trial)=(( Activations_varN_fixedW(:,trial)-(APLgains(7) )*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw)>0 ).*( Activations_varN_fixedW(:,trial)-APLgains(7)*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw);
                
            end
            
            
            for l_r=1:lrs
                
                % starting at initial random values for KC-MBon weights
                
                WopAllOdours=1*rand(n,2);
                WopAllOdoursHomog=WopAllOdours;
                WopAllOdours_Ftheta=WopAllOdours;
                WopAllOdoursHomog_Ftheta=WopAllOdours;
                
                WopAllOdoursVarW_fixedN= WopAllOdours;
                WopAllOdoursVarW_fixedN_and_theta=WopAllOdours;
                
                WopAllOdoursVarN_fixedW=WopAllOdours;
                WopAllOdoursVarN_fixedW_and_theta=WopAllOdours;
                
                WopFromPNs= 1*rand(m,2);
                %                 alpha= 0.001*(10^(l_r));
                alpha=0.000001* (10^((l_r)));
                
                c=1;
                ceq=1;
                ch=1;
                
                YHomogtemp=reshape(YHomog,n,odors,numTrials);
                YHomog_Fthetatemp=reshape(YHomog_Ftheta,n,odors,numTrials);
                Ytemp= reshape(Y,n,odors,numTrials);
                Y_Fthetatemp= reshape(Y_Ftheta,n,odors,numTrials);
                Y_varw_fixedNtemp= reshape(Y_varw_fixedN,n,odors,numTrials);
                Y_varw_fixedN_and_thetatemp= reshape(Y_varw_fixedN_and_theta,n,odors,numTrials);
                Y_varN_fixedWtemp= reshape(Y_varN_fixedW,n,odors,numTrials);
                Y_varN_fixedW_and_thetatemp=reshape(Y_varN_fixedW_and_theta,n,odors,numTrials);
                
                % rescale the KC responses to be from 0-1
                Ytemp=rescale(Ytemp);
                YHomogtemp=rescale(YHomogtemp);
                Y_Fthetatemp=rescale(Y_Fthetatemp);
                YHomog_Fthetatemp=rescale(YHomog_Fthetatemp);
                Y_varw_fixedNtemp=rescale(Y_varw_fixedNtemp);
                Y_varw_fixedN_and_thetatemp=rescale(Y_varw_fixedN_and_thetatemp);
                Y_varN_fixedWtemp=rescale(Y_varN_fixedWtemp);
                Y_varN_fixedW_and_thetatemp= rescale(Y_varN_fixedW_and_thetatemp);
                
                
                YHomogtr=YHomogtemp(:,:,1:numtrainingSamples);
                YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
                Ytr=Ytemp(:,:,1:numtrainingSamples);
                Y_Fthetatr=Y_Fthetatemp(:,:,1:numtrainingSamples);
                Y_varWtr= Y_varw_fixedNtemp(:,:,1:numtrainingSamples);
                Y_varWFtr=Y_varw_fixedN_and_thetatemp(:,:,1:numtrainingSamples);
                Y_varNtr= Y_varN_fixedWtemp(:,:,1:numtrainingSamples);
                Y_varNFtr= Y_varN_fixedW_and_thetatemp(:,:,1:numtrainingSamples);
                
                % learning the MBON- preceptron weights in the output layer in each of the 8 models: by synaptic depression
                % punishing the connection to the MBON decision (MBON avoidance & MBON approach) opposite to the odor class ('good'
                % rewarded odor and 'bad' punished odor)
                %  i.e.learning from mistakes
                
                
                for odour=1:odors
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        delta =  exp( -(alpha/mean(YHomogtr(:)))* sum(YHomogtr(:,odour,:),3) );
                        deltaWH(:,ch)= WopAllOdoursHomog(:,2).*(delta-1);
                        ch=ch+1;
                        WopAllOdoursHomog(:,2)= WopAllOdoursHomog(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(YHomogtr(:)))* sum(YHomogtr(:,odour,:),3) );
                        WopAllOdoursHomog(:,1)= WopAllOdoursHomog(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        delta =  exp( -(alpha/mean(YHomog_Fthetatr(:)))* sum(YHomog_Fthetatr(:,odour,:),3) );
                        WopAllOdoursHomog_Ftheta(:,2)= WopAllOdoursHomog_Ftheta(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(YHomog_Fthetatr(:)))* sum(YHomog_Fthetatr(:,odour,:),3) );
                        WopAllOdoursHomog_Ftheta(:,1)= WopAllOdoursHomog_Ftheta(:,1) .*delta;
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        delta =  exp( -(alpha/mean(Ytr(:)))* sum(Ytr(:,odour,:),3) );
                        deltaW(:,c)= WopAllOdours(:,2).*(delta-1);
                        c=c+1;
                        WopAllOdours(:,2)= WopAllOdours(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(Ytr(:)))* sum(Ytr(:,odour,:),3) );
                        WopAllOdours(:,1)= WopAllOdours(:,1) .*delta;
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        delta =  exp( -(alpha/mean(Y_Fthetatr(:)))* sum(Y_Fthetatr(:,odour,:),3) );
                        WopAllOdours_Ftheta(:,2)= WopAllOdours_Ftheta(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(Y_Fthetatr(:)))* sum(Y_Fthetatr(:,odour,:),3) );
                        WopAllOdours_Ftheta(:,1)= WopAllOdours_Ftheta(:,1) .*delta;
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_varWtr(:)))* sum(Y_varWtr(:,odour,:),3) );
                        WopAllOdoursVarW_fixedN(:,2)= WopAllOdoursVarW_fixedN(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(Y_varWtr(:)))* sum(Y_varWtr(:,odour,:),3) );
                        WopAllOdoursVarW_fixedN(:,1)= WopAllOdoursVarW_fixedN(:,1) .*delta;
                    end
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_varWFtr(:)))* sum(Y_varWFtr(:,odour,:),3) );
                        WopAllOdoursVarW_fixedN_and_theta(:,2)= WopAllOdoursVarW_fixedN_and_theta(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(Y_varWFtr(:)))* sum(Y_varWFtr(:,odour,:),3) );
                        WopAllOdoursVarW_fixedN_and_theta(:,1)= WopAllOdoursVarW_fixedN_and_theta(:,1) .*delta;
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_varNtr(:)))* sum(Y_varNtr(:,odour,:),3) );
                        
                        WopAllOdoursVarN_fixedW(:,2)= WopAllOdoursVarN_fixedW(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(Y_varNtr(:)))* sum(Y_varNtr(:,odour,:),3) );
                        WopAllOdoursVarN_fixedW(:,1)= WopAllOdoursVarN_fixedW(:,1) .*delta;
                    end
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_varNFtr(:)))* sum(Y_varNFtr(:,odour,:),3) );
                        WopAllOdoursVarN_fixedW_and_theta(:,2)= WopAllOdoursVarN_fixedW_and_theta(:,2) .*delta;
                        
                    else
                        delta = exp(- (alpha/mean(Y_varNFtr(:)))* sum(Y_varNFtr(:,odour,:),3) );
                        WopAllOdoursVarN_fixedW_and_theta(:,1)= WopAllOdoursVarN_fixedW_and_theta(:,1) .*delta;
                    end
                    
                    
                end
                
                %% perfromance as a function of the strictness of the decision making
                %% this strictness is dictated by C in the soft-max function.
                %% so given the same fly, same task, and after learning measure the performance as f(c)
                
                for c=1:Crange
                    
                    C=C_SoftMax*(10^c);
                    [acc,accH,acc7, acc9]=KernelTesting(C,WopAllOdours,WopAllOdoursHomog, WopAllOdoursVarW_fixedN,WopAllOdoursVarN_fixedW , classAction1,numTrials,numtrainingSamples,Ytemp,YHomogtemp, Y_varw_fixedNtemp, Y_varN_fixedWtemp);
                    
                    test_p_ra(mods,randomTrials,noiseScale, l_r,c)=acc;
                    test_p_raH(mods,randomTrials,noiseScale, l_r,c)=accH;
                    test_p_ra_varW_FixedN(mods,randomTrials,noiseScale, l_r,c)= acc7;
                    test_p_ra_varN_FixedW(mods,randomTrials,noiseScale, l_r,c)= acc9;
                    
                    [acc1,accH1,acc8,acc10]=KernelTesting_Ftheta(C,WopAllOdours_Ftheta,WopAllOdoursHomog_Ftheta,WopAllOdoursVarW_fixedN_and_theta, WopAllOdoursVarN_fixedW_and_theta,classAction1,numTrials,numtrainingSamples,Y_Fthetatemp,YHomog_Fthetatemp,Y_varw_fixedN_and_thetatemp,Y_varN_fixedW_and_thetatemp);
                    
                    test_p_ra_Fixedtheta(mods,randomTrials,noiseScale, l_r,c)=acc1;
                    test_p_raH_FixedTheta(mods,randomTrials,noiseScale, l_r,c)=accH1;
                    test_p_ra_varW_FixedN_FixedTheta(mods,randomTrials,noiseScale, l_r,c)= acc8;
                    test_p_ra_varN_FixedW_FixedTheta(mods,randomTrials,noiseScale, l_r,c)=acc10;
                    
                    
                end
                
                
                
            end
            
            
        end
        
    end
end

%cd pwd
% Models accuracies
save('test_p_ra.mat','test_p_ra');
save('test_p_raH.mat','test_p_raH');
save('test_p_ra_Fixedtheta.mat','test_p_ra_Fixedtheta');
save('test_p_raH_FixedTheta.mat','test_p_raH_FixedTheta');
save('test_p_ra_varW_FixedN_FixedTheta.mat','test_p_ra_varW_FixedN_FixedTheta')
save('test_p_ra_varW_FixedN.mat','test_p_ra_varW_FixedN')
save('test_p_ra_varN_FixedW_FixedTheta.mat','test_p_ra_varN_FixedW_FixedTheta')
save('test_p_ra_varN_FixedW.mat','test_p_ra_varN_FixedW')


% close the parallel pool
delete(gcp('nocreate'));

%------------------------------------------------END OF CODE



