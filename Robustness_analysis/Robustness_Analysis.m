%% this code is to run the robustness analysis on the main models
%% and the alternate rule for tuning the input synaptic weights (H_indiv)
%% in FigS3

clc
clear all;

% randomly 15 odor samples from the Alcohols and esters groups
% to balance the training samples in all chemical classes, 15 odor.

load('random15TerpeneAlcoholEsters.mat')

load('hallem_olsen.mat');
load('PW_given_N.mat');
load('PW_given_theta_and_n.mat');
load('W_PN_KC.mat');
% load('classAction1.mat');
load ('P_n.mat');

ll=20;
numtrainingSamples=20; %% artificial plus linearly dependent responses
lrs=7;
Crange=1;

% start parallel local thread
%p=parpool(6);
p.IdleTimeout = 10000;
parfevalOnAll(@maxNumCompThreads,0,6)

C_SoftMax=1;
NScales=1;

% each random instantiation of the network, a model fly,
% will recieve 15 odors in each chemical group, different 15 samples
% from terpenes, esters and alcohols.

for seq=1:20
    
    alcohols= AlcoOds{seq};
    esters= EstersOds{seq};
    terpenes = TerpeneOds{seq};
    
    odorsSet{seq}= [9:23,terpenes,esters,alcohols]';
    
end

tuneindexs{1}=1:15; %acids
tuneindexs{2}=16:30; %terpenes
tuneindexs{3}=31:45; %esters
tuneindexs{4}=46:60; %alcohols


for randomTrials=1:ll
    % set of odors indexes per fly: acids, 16 random esters, 16 random alcohols, terpenes
    inputsSeq_perFly= odorsSet{randomTrials};
    
    InhibitionGain= 0.0+ (0.0);
    
    APLgainP= zeros(1,2);
    
    n =2000; % number of neurons in the hidden layer
    
    m=24;  %number of dimensions in the input data
    
    clawsNo=(normrnd(6,1.7,[1 n])); %select number of claws randomly
    clawsNo(clawsNo<2)=2;
    clawsNo(clawsNo>11)=11;
    
    
    HomogClaws= ones([1,n])*6;
    
    PNsperKC = round(clawsNo.*ones(1,n));
    
    HomogPNsperKC= HomogClaws;
    
    % random KCs thresholds initializations
    ThetaMu=13;
    theta_0=abs(normrnd(ThetaMu,ThetaMu*(5.6/21.5),[n 1])); %% avoid negative values of theta
    theta_0(theta_0>70)=70;
    theta_0(theta_0<0.01)=0.01;
    
    % if this KC has only two connections, then make sure they
    % are different
    
    for i=1:n
        
        if(PNsperKC(i)==2)
            PnToKc{i} = randsample(m, PNsperKC(i), false);
            HomogPnToKc{i}= randsample(m, HomogPNsperKC(i), false);
            
        else
            PnToKc{i} = randsample(m, PNsperKC(i), true);
            HomogPnToKc{i}= randsample(m, HomogPNsperKC(i), true);
        end
        
    end
    % random initilaization of the weights
    
    thisW = zeros(m, n);
    thisW_equalizedModel_0=zeros(m,n);
    thisW_HomogModel=zeros(m,n);
    
    thisW_Kennedy_0= zeros(m,n);
    thisW_ActivityBasedComp_noxjk_0 = zeros(m,n);
    
    initial_thisW_ActivityBasedComp= zeros(m,n);
    
    
    
    
    for k=1:n
        
        for j=1:length(PnToKc{k})
            
            whichPN = PnToKc{k}(j);
            
            % pick random weight from a log normal distribution that
            % roughtly fits the Turner distribution
            
            thisWeight = exp(-0.0507+0.3527*randn(1));
            
            %while the sampled weight is zero, then
            % sample again
            while(~thisWeight)
                thisWeight = exp(-0.0507+0.3527*randn(1));
            end
            
            thisWeight_activityComp=rand(1);
            
            %while the sampled weight is zero, then
            % sample again
            while(~thisWeight_activityComp)
                thisWeight_activityComp = rand(1);
            end
            
            % sample the weights from the new fitted weights in the other script (modelling KC_PnWeights.m)
            ThetaInd= round(((theta_0(k)-0.01)/0.1)+1);
            ThetaInd(ThetaInd==0)=1;
            this_KCWeights= PW_given_theta_and_n(length(PnToKc{k})-1,ThetaInd,:);
            thisWeight_equalizedModel= randsample(W,1,'true', this_KCWeights);
            thisW(whichPN, k) = thisW(whichPN, k) + thisWeight;
            thisW_Kennedy_0(whichPN, k)= thisW_Kennedy_0(whichPN,k) + thisWeight;
            thisW_equalizedModel_0(whichPN,k)= thisW_equalizedModel_0(whichPN,k)+thisWeight_equalizedModel;
            initial_thisW_ActivityBasedComp(whichPN,k)= initial_thisW_ActivityBasedComp(whichPN,k)+ thisWeight_activityComp;
            
%             thisW_=(5+rand(1));
%             thisW_ActivityBasedComp_noxjk_0(whichPN,k)= thisW_ActivityBasedComp_noxjk_0(whichPN,k)+ thisW_;
%             
        end
    end
    
    
    for k=1:n
        
        for j=1:length(HomogPnToKc{k})
            
            whichPN_homog= HomogPnToKc{k}(j);
            
            thisWeightHomo=1; %% homogenous equal unity weights connecting KCs to PNs.
            
            thisW_HomogModel(whichPN_homog,k)= thisWeightHomo+ thisW_HomogModel(whichPN_homog,k);
            
            
        end
    end
    
    %% tune for each fly on each chemical group:
    
    for tune=1:4
        constraints=0;
        
        % tuning on tune (1:acids, 2:terpenes, 3:alchols, 4:esters)
        indexTuningOdors= tuneindexs{tune}';
        
        %train and test on others
        indexTrainingOdors= find(~ismember((1:60)',indexTuningOdors));
        
        odorsTuning=size(indexTuningOdors,1);
        odorsTuning_training= size(indexTrainingOdors,1);
        
        numTrials = 35;
        
        PN =hallem_olsen( inputsSeq_perFly(indexTuningOdors),:)';
        x=PN;
        
        PN_tune_train= hallem_olsen(inputsSeq_perFly(indexTrainingOdors),:)';
        x_tune_train= PN_tune_train;
        
        odors=size(indexTrainingOdors,1);
        
        classAction1=randsample([1:odors],round(odors/2),'false');
        
        for noiseScale=1:NScales
            
            noise=1; %((noiseScale^2)/10);
            PNtrials = zeros(24, odorsTuning, numTrials);
            PNtrials_tune_train = zeros(24, odorsTuning_training, numTrials);
            PNtrials(:,:,1) =x;
            PNtrials_tune_train(:,:,1)=x_tune_train;
            
            for t = 1:numTrials-1
                
                PNtrials(:,:,t+1) = x + ...
                    getPNStdevBhandawat(x) .* ...
                    noise.*randn(24, odorsTuning);
                
                PNtrials_tune_train(:,:,t+1) = x_tune_train + ...
                    getPNStdevBhandawat(x_tune_train) .* ...
                    noise.*randn(24, odorsTuning_training);
                
            end
            
            PNtrials(PNtrials<0)=0;
            PNtrials_tune_train(PNtrials_tune_train<0)=0;
            combinePNs = cat(2,PNtrials,PNtrials_tune_train);
            combinePNs = rescale(combinePNs,0,5);
            PNtrials = combinePNs(:,1:odorsTuning,:);
            PNtrials_tune_train = combinePNs(:,(odorsTuning+1):end,:);
%             PNtrials=rescale(PNtrials,0,5);
            
%             PNtrials_tune_train=rescale(PNtrials_tune_train,0,5);
            
            APLgainP= zeros(1,2);
            APLgainP_tune_train= zeros(1,2);
            
            tic
            C_theta=1;
            C_thetaS=1;
            C_thetaH=1;
            
            
            spmd
                % tune the P(W|n,theta) compensatory model only on the tuning odor environments
                % one chemical group
                
                if labindex==1
                    
                    thisW_equalizedModel_tune_is_train= thisW_equalizedModel_0;
                    thisW_equalizedModel=thisW_equalizedModel_0;
                    InhAbs_CLVec(1)=0;
                    InhAbsTarget=0.20;
                    InhibitionGain= 0.0+ (0.0);
                    
                    C_theta=1;
                    theta=theta_0;
                    A=zeros(n,odorsTuning*numtrainingSamples);
                    Y=zeros(n,odorsTuning*numtrainingSamples);
                    
                    ActivationsEqualizeddummy=zeros(n,odorsTuning*numtrainingSamples);
                    YEqualizeddummy=zeros(n,odorsTuning*numtrainingSamples);
                    
                    ActivationsDummy=zeros(n,odorsTuning*numtrainingSamples);
                    Activations=zeros(n,odorsTuning*numtrainingSamples);
                    Y_d=zeros(n,odorsTuning*numtrainingSamples);
                    constraints=0;
                    codingLevelDummy=[];
                    
                    while(~constraints)
                        
                        % optimization of the first error function
                        % sparsity constraint: CL(with APL blocked)~=20%
                        
                        eta=2;
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            A(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_theta.*theta) )>0 ).*( A(:,trial)-(C_theta.*theta));
                            codingLevelEqualizedDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        InhAbs_mComp=mean(codingLevelEqualizedDummy);
                        depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                        depsi1_dy(isnan(depsi1_dy))=0;
                        depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(theta,1,odorsTuning*numtrainingSamples));
                        Grad= ((InhAbs_mComp)-InhAbsTarget)*(1/(n*odorsTuning*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                        
                        C_theta=C_theta - eta.*(Grad);
                        
                        if (C_theta<0)
                            error('the scale factor in cyan model is -ve')
                        end
                        
                        %% resample the weights
                        thisW_equalizedModel=zeros(m,n);
                        for k=1:n
                            
                            for j=1:m
                                
                                if(thisW_equalizedModel_0(j,k))
                                    %% sample the weights from the new fitted weights in the other script (modelling KC_PnWeights.m)
                                    
                                    ThetaInd= round(((theta(k)-0.01)/0.1)+1);
                                    %% capping theta at 0.1..
                                    ThetaInd(ThetaInd==0)=1;
                                    this_KCWeights= PW_given_theta_and_n(PNsperKC(k)-1,ThetaInd,:);
                                    thisWeight_equalizedModel= randsample(W,1,'true', this_KCWeights);
                                    thisW_equalizedModel(j,k)= thisW_equalizedModel(j,k)+thisWeight_equalizedModel;
                                    
                                end
                                
                            end
                        end
                        
                        
                        % optimization of the second error function
                        % sparsity constraints, CL=10%
                        
                        eta_2=0.000001;
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            ActivationsEqualizeddummy(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(InhibitionGain)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)- (C_theta.*theta) )>0 ).*( ActivationsEqualizeddummy(:,trial)-InhibitionGain*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta.*theta));
                            codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                            
                        end
                        mComp=mean(codingLevelEqualizedDummy);
                        dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                        dsig_dy(isnan(dsig_dy))=0;
                        dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                        dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                        Grad_alpha= ((mComp)-0.10)*(1/(n*odorsTuning*numtrainingSamples))*(sum(dsig_dalpha(:)));
                        InhibitionGain= InhibitionGain- eta_2*(Grad_alpha);
                        
                        if (InhibitionGain<0)
                            wh=1;
                        end
                        
                        
                        % check if the constraints are satisfied
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            ActivationsDummy(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_theta.*theta))>0 ).*( ActivationsDummy(:,trial)-(C_theta.*theta));
                            codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                        end
                        
                        InhAbs_CL=mean(codingLevelDummy);
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            Activations(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            Y(:,trial)=(( Activations(:,trial)-(InhibitionGain)*repmat(sum(Activations(:,trial),1),n,1)-(C_theta.*theta))>0 ).*( Activations(:,trial)-InhibitionGain*repmat(sum(Activations(:,trial),1),n,1)-(C_theta.*theta));
                            codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        CL_=mean(codingLevelDummy);
                        constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.1) &( abs(CL_-0.10)<0.01 );
                        
                        InhAbs_CLVec(end+1)=InhAbs_CL;
                        
                    end
                    CLevelP(1)=CL_;
                    INHAbs_CLP(1)=InhAbs_CL;
                    theta=(C_theta.*theta);
                    InhibitionGain_tune_is_train= 0.0+ (0.0);
                    
                    % tuning the same network bUT using the same odors of
                    % training, a familiar environment.
                    
                    theta_tune_is_train= theta_0;
                    C_theta_1=1;
                    
                    [thisW_equalizedModel_tune_is_train,theta_tune_is_train, InhibitionGain_tune_is_train,CL_, InhAbs_CL] = Tune_Pw_G_N_T_Model_on_trainingData(PNsperKC,PW_given_theta_and_n,W,thisW_equalizedModel_tune_is_train,...
                        PNtrials_tune_train, theta_tune_is_train,InhibitionGain_tune_is_train,odors,numtrainingSamples,C_theta_1 );
                    
                    CLevelP_tune_train(1)=CL_;
                    INHAbs_CLP_tune_train(1)=InhAbs_CL;
                    %APLgainP_tune_train(2)=APLgainP_tune_train_H;
                    
                    InhibitionGain_tune_train(1)=InhibitionGain_tune_is_train;
                    
                    
                    theta_tune_train{1}=theta_tune_is_train;
                    thisW_equalizedModel_tune_train{1} = thisW_equalizedModel_tune_is_train;
                    
                    
                    %% for the random model:
                    
                elseif labindex==2
                    
                    T=10;
                    thetaS_0=abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
                    thetaS_0(thetaS_0>70)=70;
                    thetaS_0(thetaS_0<0.01)=0.01;
                    constraints=0; %% again for the variable theta, random and homog models
                    A=zeros(n,odorsTuning*numtrainingSamples);
                    Y=zeros(n,odorsTuning*numtrainingSamples);
                    ActivationsEqualizeddummy=zeros(n,odorsTuning*numtrainingSamples);
                    YEqualizeddummy=zeros(n,odorsTuning*numtrainingSamples);
                    Activations_S=zeros(n,odorsTuning*numtrainingSamples);
                    Activation=zeros(n,odorsTuning*numtrainingSamples);
                    Y_S=zeros(n,odorsTuning*numtrainingSamples);
                    
                    eta=1;
                    
                    thetaS=thetaS_0;
                    
                    while(~constraints)
                        
                        % first sparsity constraint, CL(with APL
                        % blocked)~=20%
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            A(:,trial) = thisW'*PNtrials(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_thetaS .*thetaS))>0 ).*( A(:,trial)-(C_thetaS .*thetaS));
                            codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        InhAbs_mSimp=mean(codingLevelDummy);
                        depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                        depsi1_dy(isnan(depsi1_dy))=0;
                        depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(thetaS,1,odorsTuning*numtrainingSamples));
                        Grad= ((InhAbs_mSimp)-0.20)*(1/(n*odorsTuning*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                        C_thetaS=C_thetaS - eta.*(Grad);
                        
                        if (C_thetaS<0)
                            error('the scale factor in the random model is -ve')
                        end
                        
                        % second sparsity constraint,CL~=10%
                        eta_2=0.0000001;
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            ActivationsEqualizeddummy(:,trial) = thisW'*PNtrials(:,trial);
                            YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(1))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaS .*thetaS))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(1)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaS .*thetaS));
                            codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                            
                        end
                        CLRand=mean(codingLevelEqualizedDummy);
                        dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                        dsig_dy(isnan(dsig_dy))=0;
                        dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                        dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                        Grad_alpha= ((CLRand)-0.10)*(1/(n*odorsTuning*numtrainingSamples))*(sum(dsig_dalpha(:)));
                        APLgainP(1)= APLgainP(1)- eta_2*(Grad_alpha);
                        
                        % check if constraints are satisfied
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            Activations_S(:,trial) = thisW'*PNtrials(:,trial);
                            Y_S(:,trial)=(( Activations_S(:,trial)-(C_thetaS.*thetaS) )>0 ).*( Activations_S(:,trial)-(C_thetaS.*thetaS));
                            codingLevelDummy(trial)=  (sum(Y_S(:,trial)>0,1)/n);
                        end
                        InhAbs_CL=mean(codingLevelDummy);
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            Activation(:,trial) = thisW'*PNtrials(:,trial);
                            Y(:,trial)=(( Activation(:,trial)-(APLgainP(1))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS))>0 ).*( Activation(:,trial)-APLgainP(1)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS));
                            codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        CL_=mean(codingLevelDummy);
                        
                        constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.1 ) &( abs(CL_-0.10)<0.01 );
                        
                        
                    end
                    CLevelP(1)=CL_;
                    INHAbs_CLP(1)=InhAbs_CL;
                    thetaS= (C_thetaS.*thetaS);
                    
                    C_thetaS_1=1;
                    thetaS_tune_train=thetaS_0;
                    
                    [thetaS_tune_train, APLgainP_tune_train_S,CL_, InhAbs_CL] = Tune_Random_Model_on_trainingData(thisW,...
                        PNtrials_tune_train, thetaS_tune_train, APLgainP_tune_train(1),odorsTuning_training,numtrainingSamples,C_thetaS_1);
                    
                    CLevelP_tune_train(1)=CL_;
                    INHAbs_CLP_tune_train(1)=InhAbs_CL;
                    APLgainP_tune_train(1)=APLgainP_tune_train_S;
                    
                    
                    
                    
                elseif labindex==3
                    
                    % the Homogenous model
                    
                    constraints=0;
                    thetaH_Ftheta_0=10+rand(1);
                    A=zeros(n,odorsTuning*numtrainingSamples);
                    Y=zeros(n,odorsTuning*numtrainingSamples);
                    ActivationsEqualizeddummy=zeros(n,odorsTuning*numtrainingSamples);
                    YEqualizeddummy=zeros(n,odorsTuning*numtrainingSamples);
                    Activations_H=zeros(n,odorsTuning*numtrainingSamples);
                    Y_H=zeros(n,odorsTuning*numtrainingSamples);
                    Activation=zeros(n,odorsTuning*numtrainingSamples);
                    Y=zeros(n,odorsTuning*numtrainingSamples);
                    HomoFtheta_Inh=0;
                    thetaH_Ftheta= thetaH_Ftheta_0;
                    
                    while(~constraints)
                        
                        eta=1;
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            A(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_thetaH.*thetaH_Ftheta))>0 ).*( A(:,trial)-(C_thetaH.*thetaH_Ftheta));
                            codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        InhAbs_mHomog=mean(codingLevelDummy);
                        depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                        depsi1_dy(isnan(depsi1_dy))=0;
                        depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaH_Ftheta,n,odorsTuning*numtrainingSamples));
                        Grad= ((InhAbs_mHomog)-0.20)*(1/(n*odorsTuning*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                        C_thetaH=C_thetaH - eta.*(Grad);
                        
                        if (C_thetaH<0)
                            error('the scale factor in black model is -ve')
                        end
                        
                        eta_2=0.0000001;
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            ActivationsEqualizeddummy(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                            YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(APLgainP(2))*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaH.*thetaH_Ftheta))>0 ).*( ActivationsEqualizeddummy(:,trial)-APLgainP(2)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_thetaH.*thetaH_Ftheta));
                            codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                            
                        end
                        CLAllFixed=mean(codingLevelEqualizedDummy);
                        dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                        dsig_dy(isnan(dsig_dy))=0;
                        dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                        dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                        Grad_alpha= ((CLAllFixed)-0.10)*(1/(n*odorsTuning*numtrainingSamples))*(sum(dsig_dalpha(:)));
                        APLgainP(2)= APLgainP(2)- eta_2*(Grad_alpha);
                        
                        % check if the constraints are satisfied
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            Activations_H(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                            Y_H(:,trial)=(( Activations_H(:,trial)-(C_thetaH.*thetaH_Ftheta))>0 ).*( Activations_H(:,trial)-(C_thetaH.*thetaH_Ftheta));
                            codingLevelDummy(trial)=  (sum(Y_H(:,trial)>0,1)/n);
                        end
                        InhAbs_CL=mean(codingLevelDummy);
                        
                        for trial = 1:(odorsTuning*numtrainingSamples)
                            
                            Activation(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                            Y(:,trial)=(( Activation(:,trial)-(APLgainP(2))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH.*thetaH_Ftheta))>0 ).*( Activation(:,trial)-APLgainP(2)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH.*thetaH_Ftheta));
                            codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        CL_=mean(codingLevelDummy);
                        
                        constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.1) &( abs(CL_-0.10)<0.01 );
                        
                        
                    end
                    CLevelP(1)=CL_;
                    INHAbs_CLP(1)=InhAbs_CL;
                    thetaH_Ftheta= (C_thetaH.*thetaH_Ftheta);
                    
                    
                    
                    thetaH_Ftheta_tune_train= thetaH_Ftheta_0;
                    C_thetaH_1=1;
                    
                    [thetaH_Ftheta_tune_train, APLgainP_tune_train_H,CL_, InhAbs_CL] = Tune_Homogenous_Model_on_trainingData(thisW_HomogModel,...
                        PNtrials_tune_train, thetaH_Ftheta_tune_train, APLgainP_tune_train(2),odorsTuning_training,numtrainingSamples, C_thetaH_1  );
                    
                    CLevelP_tune_train(1)=CL_;
                    INHAbs_CLP_tune_train(1)=InhAbs_CL;
                    APLgainP_tune_train(2)=APLgainP_tune_train_H;
                    
                    
                end
            end
            Clevels=[CLevelP{1},CLevelP{2},CLevelP{3}];
            INHAbs_CL=[INHAbs_CLP{1},INHAbs_CLP{2},INHAbs_CLP{3}];
            
            Clevels_tune_is_train=[CLevelP_tune_train{1},CLevelP_tune_train{2},CLevelP_tune_train{3}];
            INHAbs_CL_tune_is_train=[INHAbs_CLP_tune_train{1},INHAbs_CLP_tune_train{2},INHAbs_CLP_tune_train{3}];
            
            
            for i=1:2
                temp=APLgainP{i+1};
                APLgains(i)= temp(i);
                
                temp1=APLgainP_tune_train{i+1};
                APLgains_tune_is_train(i)= temp1(i);
                
            end
            
            InhibitionGain=InhibitionGain{1};
            InhibitionGain_tune_is_train=InhibitionGain_tune_train{1};
            theta=theta{1};
            thetaH_Ftheta=thetaH_Ftheta{3};
            thetaS=thetaS{2};
            thisW_equalizedModel=thisW_equalizedModel{1};
            
            theta_tune_is_train=theta_tune_train{1};
            if iscell(theta_tune_is_train)
                if length(theta_tune_is_train)==1
                    theta_tune_is_train = theta_tune_is_train{1};
                else
                    error('something is wrong with theta_tune_is_train')
                end
            end
            
            thetaH_Ftheta_tune_is_train=thetaH_Ftheta_tune_train{3};
            thetaS_tune_is_train=thetaS_tune_train{2};
            thisW_equalizedModel_tune_is_train=thisW_equalizedModel_tune_train{1};
            if iscell(thisW_equalizedModel_tune_is_train)
                if length(thisW_equalizedModel_tune_is_train)==1
                    thisW_equalizedModel_tune_is_train = thisW_equalizedModel_tune_is_train{1};
                else
                    error('something is wrong with thisW_equalizedModel_tune_is_train')
                end
            end
            
            t_cyan_red_black=toc;
            
            tic
            %% optimization for the activity dependent models
            
            %% tuning the input PN to KC synaptic weights as in our
            %% main model, i.e. all PN input weights to a KC are scaled homogenously
            
            T=5;
            theta_comp2_noxjk_0= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_comp2_noxjk_0(theta_comp2_noxjk_0>70)=70;
            theta_comp2_noxjk_0(theta_comp2_noxjk_0<0.01)=0.01;
            
            theta_comp2_noxjk=theta_comp2_noxjk_0;
            theta_comp2_noxjk_tuneis_train=theta_comp2_noxjk_0;
            thisW_ActivityBasedComp_noxjk=initial_thisW_ActivityBasedComp;
            thisW_ActivityBasedComp_noxjk_tuneis_train=initial_thisW_ActivityBasedComp;
            Activations =zeros(n,odorsTuning*numtrainingSamples);
            ActivationsDummy= zeros(n, odorsTuning*numtrainingSamples);
            Conn=zeros(m,n);
            Conn(find(thisW_ActivityBasedComp_noxjk))=1;
            mask= zeros (m,n);
            mask(find(thisW_ActivityBasedComp_noxjk))=1;
            C_=1;
            APLgains_noxjk=0;
            A0=(0.51)*ones(n,1);
            epsilon= A0(1)*0.06;
            Inhabs_CLV=[];
            CLV=[];
            conditions=0; %%
            A=zeros(n,odorsTuning*numtrainingSamples);
            Y_d=zeros(n,odorsTuning*numtrainingSamples);
            Y_=[];
            codingLevelDummy=[];
            iter_till_exit=0;
            t=1;
            while( (~conditions)&&(iter_till_exit<10000) )
                
                iter_till_exit=iter_till_exit+1;
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    A(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_.*theta_comp2_noxjk) )>0 ).*( A(:,trial)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_comp2_noxjk,1,odorsTuning*numtrainingSamples));
                eta=10;
                Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
                C_=C_ - (eta*Grad);
                if (C_<0)
                    error('the scale factor in the blue model is -ve!!')
                end
                
                eta_2=0.00000005;
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk))>0 ).*( Activations(:,trial)-APLgains_noxjk*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(Activations,1);
                
                if ( any(isinf(dAct_dalpha)) )
                    stopp=1;
                end
                dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_)-0.10)*(1/(n*odorsTuning*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains_noxjk= APLgains_noxjk- eta_2*(Grad_alpha);
                
                
                %% tuning the KC input weights to reach the average target activity A0
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial );
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk))>0 ).*( Activations(:,trial)-APLgains_noxjk*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk));
                end
                avgAKcs=mean(Y_,2);
                errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);
                thisW_ActivityBasedComp_noxjk= thisW_ActivityBasedComp_noxjk-(0.05).*((1.*(mask.*errorInActivity)));
                
                %catch the -ve weights values
                if (~isempty(find(isinf(thisW_ActivityBasedComp_noxjk) )))
                    g=1;
                end
                thisW_ActivityBasedComp_noxjk(find(thisW_ActivityBasedComp_noxjk<0))=0;
                
                % check if the constraints are satisfied
                for trial = 1:(odorsTuning*numtrainingSamples)
                    ActivationsDummy(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_comp2_noxjk))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    Activations(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk))>0 ).*( Activations(:,trial)-APLgains_noxjk*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                if mod(t,10)==0
                    disp('blue');
                    disp(CL_)
                    disp(nnz(abs(avgAKcs-A0)<epsilon))
                end
                t=t+1;
                TunedKCs= size(find((abs(avgAKcs-A0)<epsilon)));
                conditions= TunedKCs(1)>=1995 &( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01 );
                
            end
            
            % rescue the initial states of the weights and fill the
            % parameters to be tuned with NANs, failed to converge
             if(iter_till_exit<10000)
                 theta_comp2_noxjk=(C_.*theta_comp2_noxjk); 
                 toc
                 disp('blue model done');
                
             else
                 disp('blue model failed to converge on unseen data');
                 thisW_ActivityBasedComp_noxjk(find(thisW_ActivityBasedComp_noxjk))=nan;
                 theta_comp2_noxjk=ones(2000,1).*nan;
                 APLgains_noxjk=nan;
                 save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_blue_unseenData'),'theta_comp2_noxjk_0','initial_thisW_ActivityBasedComp');
                 
             end

            % rescue the initial states of the weights and fill the
            % parameters to be tuned with NANs, failed to converge
          
            %% tuning the input weights on the same odors that would
            % be in the training phase; a familiar environment
           
            Tune_darkBlueModel_on_trainingData;
            
            % rescue the initial states of the weights and fill the
            % parameters to be tuned with NANs, failed to converge
             if(iter_till_exit<10000)
                theta_comp2_noxjk_tuneis_train=(C_.*theta_comp2_noxjk_tuneis_train);
                toc
                disp('finished tune noxjk on training');
             else
                 disp('blue model failed to converge on training data');
                 thisW_ActivityBasedComp_noxjk_tuneis_train(find(thisW_ActivityBasedComp_noxjk_tuneis_train))=nan;
                 theta_comp2_noxjk_tuneis_train=ones(2000,1).*nan;
                 APLgains_noxjk_tuneis_train=nan;
                 save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_blue_trainingData'),'theta_comp2_noxjk_0','initial_thisW_ActivityBasedComp');
                 
             end

            % rescue the initial states of the weights and fill the
            % parameters to be tuned with NANs, failed to converge
            
          %% tuning the input weights on the same odors that would
            % be in the training phase; a familiar environment
            
            tic
            %% tuning the input PN to KC synaptic weights using the rule
            %% with the extra <xjk>_k factor, i.e. each PN-KC input is tuned individually
            
            T=5;
            theta_comp2_0= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_comp2_0(theta_comp2_0>70)=70;
            theta_comp2_0(theta_comp2_0<0.01)=0.01;
            
            thisW_ActivityBasedComp= initial_thisW_ActivityBasedComp;
            thisW_ActivityBasedComp_tune_is_train= initial_thisW_ActivityBasedComp;
            
            APLgains(3)=0;
            APLgains_tune_is_train(3)=0;
            Activations =zeros(n,odorsTuning*numtrainingSamples);
            ActivationsDummy= zeros(n, odorsTuning*numtrainingSamples);
            A=zeros(n,odorsTuning*numtrainingSamples);
            Y_d=zeros(n,odorsTuning*numtrainingSamples);
            Y_= zeros(n,odorsTuning*numtrainingSamples);
            codingLevelDummy=[];
            
            % target average activity
            A0=(0.51)*ones(n,1);
            % error or premissible tolerance from the target
            % average activity A0
            epsilon= A0(1)*0.07;
            
            theta_comp2= theta_comp2_0;
            conditions=0;
            Conn= zeros (m,n);
            Conn(find(thisW_ActivityBasedComp))=1;
            C_=1;
            iter_till_exit=0;
            t=1;
            while(~conditions && (iter_till_exit<10000))
                
                iter_till_exit=iter_till_exit+1;
                eta=10;
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    A(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_.*theta_comp2))>0 ).*( A(:,trial)-(C_.*theta_comp2));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevelDummy);
                depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_comp2,1,odorsTuning*numtrainingSamples));
                Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
                C_=C_ - eta*(Grad);
                
                if (C_<0)
                    error('the scale factor in the blue model is -ve!!')
                end
                
                eta_2=0.00000005;
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(3))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2))>0 ).*( Activations(:,trial)-APLgains(3)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(Activations,1);
                
                if ( any(isinf(dAct_dalpha)) )
                    stopp=1;
                end
                dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_)-0.10)*(1/(n*odorsTuning*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(3)= APLgains(3)- eta_2*(Grad_alpha);
                
                
                %% tuning the input weights for activity equalization
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial );
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(3))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2))>0 ).*( Activations(:,trial)-APLgains(3)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2));
                end
                avgAKcs=mean(Y_,2);
                errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);
                Conn2=repmat(Conn,1,1,odorsTuning*numtrainingSamples);
                Xik=reshape(PNtrials(:,1:odorsTuning*numtrainingSamples),m,1,odorsTuning*numtrainingSamples);
                Xik_=repmat(Xik,1,n,1);
                filt_=  ( ((1-(InhibitionGain/n)) .* (Xik_.* Conn2)) ) ;
                filt_Xjm= mean(filt_,3);
                mask=filt_Xjm;
                
                thisW_ActivityBasedComp= thisW_ActivityBasedComp-(0.15).*((1.*(mask.*errorInActivity)));
                
                %catch the -ve weights values
                if (~isempty(find(isinf(thisW_ActivityBasedComp) )))
                    g=1;
                end
                thisW_ActivityBasedComp(find(thisW_ActivityBasedComp<0))=0;
                
                % constarints satisfaction check:
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    ActivationsDummy(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_comp2))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_comp2));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(3))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2))>0 ).*( Activations(:,trial)-APLgains(3)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                avgAKcs= mean(Y_,2);
                avgact_trace(:,t) = avgAKcs;
                if mod(t,10)==0
                    disp('dark blue');
                    disp(t)
                    disp(CL_)
                    disp(InhAbs_CL);
                    disp(max(avgAKcs))
                    disp(min(avgAKcs))
                    disp(nnz(abs(avgAKcs-A0)<epsilon))
                end
                t=t+1;
                TunedKCs= size(find((abs(avgAKcs-A0)<epsilon)));
                conditions= TunedKCs(1)>=1995 &( abs( (InhAbs_CL/CL_) - 2.0)<0.1 ) &( (abs(CL_-0.10)) <=0.01 );
                
            end
            Clevels(end+1)=CL_;
            INHAbs_CL(end+1)=InhAbs_CL;
            
             % rescue the initial states of the weights and fill the
            % parameters to be tuned with NANs, failed to converge
             if(iter_till_exit<10000)
                 theta_comp2=(C_.*theta_comp2);                
                 toc
                 disp('dark blue model done');
                 
             else
                 disp('dark blue model failed to converge on unseen data');
                 thisW_ActivityBasedComp(find(thisW_ActivityBasedComp))=nan;
                 theta_comp2=ones(2000,1).*nan;
                 APLgains(3)=nan;
                 save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_darkBlue_unseenData'),'theta_comp2_0','initial_thisW_ActivityBasedComp');
                 
             end

            % rescue the initial states of the weights and fill the
            % parameters to be tuned with NANs, failed to converge
          
            %% tuning the input weights on the same odors that would
            % be in the training phase; a familiar environment
            C_1=1;
            theta_comp2_tune_is_train= theta_comp2_0;
            tic
            
            Tune_Homeostatic_Model_on_trainingData;
            
            if(iter_till_exit<10000)
                
            theta_comp2_tune_is_train=(C_1.*theta_comp2_tune_is_train);
            toc
            disp('finished Homeostatic model on training')
                
            else
                disp('dark blue model failed to converge on training data');
                 thisW_ActivityBasedComp_tune_is_train(find(thisW_ActivityBasedComp_tune_is_train))=nan;
                 theta_comp2_tune_is_train=ones(2000,1).*nan;
                 APLgains_tune_is_train(3)=nan;
                 save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_darkBlue_trainingData'),'theta_comp2_0','initial_thisW_ActivityBasedComp');
                 
                
            end
            
            % sanity check
            if (thetaH_Ftheta<0)
                error('the theta in the homogeneous model cant be negative')
            end
            
            
            tic
            %% theta tuned for activity equalization
            
            Activations =zeros(n,odorsTuning*numtrainingSamples);
            ActivationsDummy= zeros(n, odorsTuning*numtrainingSamples);
            
            % pre-scaling the PN-KC weight matrix to achieve the
            % target activity level.
            if(tune~=3)
                thisW_Kennedy=5.5.*thisW_Kennedy_0;
            else
                thisW_Kennedy= 4.0.*thisW_Kennedy_0;
            end
            
            APLgains(4)=0;
            APLgains_tune_is_train(4)=0;
            
            if tune==3
                theta_Activity_homeo_0=  10+ 1.*rand(2000,1);
            else
                theta_Activity_homeo_0=  5+ 1.*rand(2000,1); %% avoid negative values of theta
            end
            conditions=0;
            theta_Activity_homeo=theta_Activity_homeo_0;
            A=zeros(n,odorsTuning*numtrainingSamples);
            Y_d=zeros(n,odorsTuning*numtrainingSamples);
            Y_= zeros(n,odorsTuning*numtrainingSamples);
            if tune==3
                A0=0.26.*ones(n,1);
            else
                A0=(0.3).*ones(n,1);
            end
            epsilon= A0(1)*0.07;
            C_=1;
            t=1;
            iterr=1;
            eta_gradAct_theta_0=0.15;%0.05;
            drop=0.7;
            iterDrop=1000;
            cw=1;
            
            if tune==3
                eta_0 = 0.01;
            else
                eta_0=0.05;
            end
            drop_1=0.5;
            iterrDrop_1=1000;
            
            codingLevelDummy=[];
            iter_till_exit=0;
            
            while(~conditions &&(iter_till_exit<10000))
                
                iter_till_exit=iter_till_exit+1;
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    A(:,trial) = (cw.*thisW_Kennedy)'*PNtrials(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_.*theta_Activity_homeo))>0 ).*( A(:,trial)-(C_.*theta_Activity_homeo));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevelDummy);
                depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_Activity_homeo,1,odorsTuning*numtrainingSamples));
                eta= eta_0; %*(drop_1^(floor(iterr/iterrDrop_1)));
                Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
                C_=C_- (eta*Grad);
                if (C_<0)
                    error('the scale factor in our magenta model is -ve')
                end
                
                %% CL=10%, sparsity constraint
                % decaying learning rate: as more iterations, it means
                % we are closer to the minima and hence should
                % update the APL_gain less.
                
                eta_2=0.000000005; %*(0.7^(floor(iterr/1000)));
                for trial = 1:(odorsTuning*numtrainingSamples)
                    Activations(:,trial) = (cw.*thisW_Kennedy)'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(4))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_Activity_homeo))>0 ).*( Activations(:,trial)-APLgains(4)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_Activity_homeo));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevelDummy);
                dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(Activations,1);
                
                if ( any(isinf(dAct_dalpha)) )
                    
                    stopp=1;
                end
                
                dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_)-0.10)*(1/(n*odorsTuning*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(4)= APLgains(4)- eta_2*(Grad_alpha);
                
                if(APLgains(4)<0)
                    see=1;
                end
                
                %% now tune theta for equalizing average activity level.
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = (cw.*thisW_Kennedy)'*PNtrials(:,trial );
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(4))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_Activity_homeo))>0 ).*( Activations(:,trial)- APLgains(4)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_Activity_homeo));
                    
                    %can yield +ve and -ve values
                    Act_minus_inhibition(:,trial)=(Activations(:,trial)-(APLgains(4))*repmat(sum(Activations(:,trial),1),n,1));
                    
                end
                eta_gradAct_theta=eta_gradAct_theta_0; %*(drop^floor(iterr/iterDrop));
                avgAKcs=mean(Y_,2);
                
                errorInActivity=(avgAKcs-A0);
                
                theta_Activity_homeo= theta_Activity_homeo-(((eta_gradAct_theta).*((-1.*C_.*(errorInActivity)))));
                theta_Activity_homeo(theta_Activity_homeo<0)=0;
                
                %% check the constraints
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    ActivationsDummy(:,trial) = (cw.*thisW_Kennedy)'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)- (C_.*theta_Activity_homeo))>0 ).*( ActivationsDummy(:,trial)- (C_.*theta_Activity_homeo));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = (cw.*thisW_Kennedy)'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(4))*repmat(sum(Activations(:,trial),1),n,1)- (C_.*theta_Activity_homeo))>0 ).*( Activations(:,trial)-APLgains(4)*repmat(sum(Activations(:,trial),1),n,1)- (C_.*theta_Activity_homeo));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                avgAKcs=mean(Y_,2);
                
                
                avgact_trace(:,iterr)=avgAKcs;
                t=t+1;
                if(any(theta_Activity_homeo==0))
                    key=1;
                end
                
                %specifically ONLY for this model:
                % I can accept a solution with at least 1995 KCs (well tuned) or more
                % allowing straying of 5KCs away from the target
                % average activity
                TunedKCs=size(find(abs(avgAKcs-A0)<epsilon));
                
                conditions= TunedKCs(1)>=1995 &( round(abs( round((InhAbs_CL/CL_),1) - 2.0),1) <=0.3 ) &( (abs(CL_-0.10)) <=0.015 );
                if mod(iterr,10)==0
                    disp('magenta');
                    %disp(eta)
                    %disp(eta_2)
                    disp(iterr)
                    disp(CL_)
                    disp(InhAbs_CL)
                    disp( TunedKCs(1))
                    
                end
                
                iterr=iterr+1;
                
            end
            
             if(iter_till_exit<10000)
                 
                disp('magenta model on unseen data done')
                theta_Activity_homeo=C_.*theta_Activity_homeo;
                thisW_Kennedy=cw.*thisW_Kennedy;
            
            else
                disp('magenta model failed to converge on unseen data');
                 theta_Activity_homeo=ones(2000,1).*nan;
                 APLgains(4)=nan;
                 save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_magenta_unseenData'),'thisW_Kennedy','theta_Activity_homeo_0');
                 
                
            end
            
            % tuning the same model on the familiar environment;
            % the odors that will be then used for training too.
            
            C_1=1;
            if tune==3
                theta_Activity_homeo_0= 5+ 1.*rand(2000,1);
            else
            theta_Activity_homeo_0= 2+ 1.*rand(2000,1);
            end
            theta_Activity_homeo_tune_is_train=theta_Activity_homeo_0;
           
            Tune_Theta_Homeo_Model_on_trainingData;
            
            if(iter_till_exit<10000)
                theta_Activity_homeo_tune_is_train=C_1.*theta_Activity_homeo_tune_is_train;
                toc
                disp('finished tune thetahomeo on training');

            else

                 disp('magenta model failed to converge on training data');
                 theta_Activity_homeo_tune_is_train=ones(2000,1).*nan;
                 APLgains_tune_is_train(4)=nan;
                 save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_magenta_trainingData'),'thisW_Kennedy','theta_Activity_homeo_0');
            end
           
            
            
            %% tuning the input inhibitory weights, APL->KC input weights
            tic
            APLtrajectory=zeros(2000,1);
            CLtrajectory=[];
            AvgAKCstrajectory=zeros(2000,1);
            C_=10;
            T=1;
            theta_inhibitionPlast_0= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_inhibitionPlast_0(theta_inhibitionPlast_0>70)=70;
            theta_inhibitionPlast_0(theta_inhibitionPlast_0<0.01)=0.01;
            
            if (tune~=3)
                thisW_ActivityBasedComp_inhibitionPlast= 5.*thisW;
            else
                thisW_ActivityBasedComp_inhibitionPlast= 5.*thisW;
            end
            
            Activations =zeros(n,odorsTuning*numtrainingSamples);
            ActivationsDummy= zeros(n, odorsTuning*numtrainingSamples);
            
            APLgains_model6=zeros(2000,1);
            APLgains_model6_tune_is_train=zeros(2000,1);
            
            iterr=1;
            
            A0=(0.32).*ones(n,1);
            
            epsilon= A0(1)*0.07;
            
            conditions=0;
            
            A=zeros(n,odorsTuning*numtrainingSamples);
            Y_d=zeros(n,odorsTuning*numtrainingSamples);
            Y_= zeros(n,odorsTuning*numtrainingSamples);
            codingLevelDummy=[];
            
            theta_inhibitionPlast=theta_inhibitionPlast_0;
            iter_till_exit=0;
            t=1;
            
            while(~conditions &&(iter_till_exit<10000) )
                iter_till_exit=iter_till_exit+1;
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    A(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_.*theta_inhibitionPlast))>0 ).*( A(:,trial)-(C_.*theta_inhibitionPlast));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevelDummy);
                
                depsi1_dy=(exp(0.91.*Y_d)./((1+exp(0.91.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_inhibitionPlast,1,odorsTuning*numtrainingSamples));
                
                if tune==3
                    eta=50;
                else
                    eta=15; %5; %*(0.7^(floor(iterr/1000))); % decaying learning rate
                end
                
                Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
                
                C_=C_ -(eta.* mean(Y_d(:)).*Grad);
                
                if (C_<0)
                    error('the scale factor in the green model is -ve!')
                    
                end
                
                % tuning inhibition weights for both: (a)
                % network sparsity level=10% (b) equalizing
                % average activty levels for all KCs
                % (homoestasis-like
                % aspect of the model)
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    Activations(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_model6.*(sum(Activations(:,trial),1)))-(C_.*theta_inhibitionPlast))>0 ).*( Activations(:,trial)-(APLgains_model6.*(sum(Activations(:,trial),1)))-(C_.*theta_inhibitionPlast));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                avgAKcs=mean(Y_,2);
                errorInActivity=avgAKcs-A0;
                dsig_dy=(exp(0.91.*Y_)./((1+exp(0.91.*Y_)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                xk_=repmat(sum(Activations),n,1);
                D_yj_alphaj= (-1.*mean(xk_,2));
                dAct_dalpha= sum(Activations,1);
                
                dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                
                dYik_dalphai= -(repmat(dAct_dalpha,n,1));
                
                if tune==3
                    eta_01 = 0.002;
                else
                    eta_o1= 0.01;
                end
                eta_o2= 5e-7;%1e-6; %5e-7;
                
                Grad_alpha1= ( ((eta_o1)) .*((CL_)-0.10)*(1/(n*odorsTuning*numtrainingSamples)).*(sum(dsig_dalpha,2)) ) ;
                Grad_alpha2=( eta_o2.* (errorInActivity).*(D_yj_alphaj) ) ;
                Grad_alpha= Grad_alpha1+Grad_alpha2 ;
                APLgains_model6= APLgains_model6- 0.001.*((Grad_alpha));
                
                
                %% check constraints
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    ActivationsDummy(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_inhibitionPlast))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_inhibitionPlast));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odorsTuning*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_model6.*sum(Activations(:,trial),1))-(C_.*theta_inhibitionPlast))>0 ).*( Activations(:,trial)-(APLgains_model6.*sum(Activations(:,trial),1))-(C_.*theta_inhibitionPlast));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                avgAKcs=mean(Y_,2);
                
                TunedKCs=size(find(abs(avgAKcs-A0)<epsilon));
                conditions= TunedKCs(1)>=1995  &( (abs(round(CL_,3)-0.10)) <=0.015) & ( round( abs( ((InhAbs_CL/CL_)) - 2.0),1) <=0.2 );
                
                if mod(iterr,10)==0
                    disp('green');
                    disp(iterr)
                    disp(CL_)
                    disp(InhAbs_CL)
                    disp(nnz(abs(avgAKcs-A0)<epsilon))
                end
                
                iterr=iterr+1;
                avgact_trace(:,t)=avgAKcs;
                t=t+1;
                
            end
            Clevels(end+1)=CL_;
            INHAbs_CL(end+1)=InhAbs_CL;
            
            if (iter_till_exit<10000)
                disp('finished green model on unseen data')
                theta_inhibitionPlast=(C_.*theta_inhibitionPlast);
                toc
            
            else
                disp('failed to converge on green model unseen data')
                theta_inhibitionPlast=ones(2000,1).*nan;
                APLgains_model6=ones(2000,1).*nan;
                save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_green_unseenData'),'thisW_ActivityBasedComp_inhibitionPlast','theta_inhibitionPlast_0');

            end
            
            % tuning the same model on the familiar environment;
            % the odors that will be then used for training too.
            C_1=1;
            theta_inhibitionPlast_tune_is_train=theta_inhibitionPlast_0;
            Tune_InhPlast_Model_on_trainingData;
            
            if (iter_till_exit<10000)
                disp('finished green model on training data')
               Clevels_tune_is_train (end+1)=CL_;
               INHAbs_CL_tune_is_train(end+1)=InhAbs_CL;
               theta_inhibitionPlast_tune_is_train=(C_1.*theta_inhibitionPlast_tune_is_train);
                toc

            else
                disp('failed to converge on green model training data')
                theta_inhibitionPlast_tune_is_train=ones(2000,1).*nan;
                APLgains_model6_tune_is_train=ones(2000,1).*nan;
                save(strcat('fly',num2str(randomTrials),'_failed_to_converge_on_tune',num2str(tune),'_model_green_trainingData'),'thisW_ActivityBasedComp_inhibitionPlast','theta_inhibitionPlast_0');

            end
            
            
            
            % if the loops broke without the constraints being statisfied. then
            % print an error message
            
            if( (any(APLgains<0)) || (any(APLgains_tune_is_train<0))  )
                
                error( 'constraints cant be statisfied ' )
                
            end
            
            %% uncomment if to: save the tuned parameters for each fly, in each tuning phase,
            %% for all the models:
            
            if(odors==size(indexTrainingOdors,1))
                save( strcat('robustnessTest_tunedSubgroup_trainedAllOtherGroups_BalancedEsters_and_Alcohols_tuneNo_',num2str(tune),[' _fly_wNoise',num2str(randomTrials),num2str(noiseScale)]) , 'thisW_HomogModel','APLgains',...
                    'thetaH_Ftheta','thisW_Kennedy', 'thisW','theta', 'thetaS', 'theta_Activity_homeo','InhibitionGain','thisW_equalizedModel','PNtrials_tune_train','PNtrials','theta_comp2','thisW_ActivityBasedComp','theta_inhibitionPlast','APLgains_model6','thisW_ActivityBasedComp_inhibitionPlast','theta_comp2_noxjk','thisW_ActivityBasedComp_noxjk','APLgains_noxjk');
                
                save( strcat('robustnessTest_tuned_and_trained_onSameGroups_BalancedEsters_and_Alcohols_tuneNo_',num2str(tune),[' _fly_wNoise',num2str(randomTrials),num2str(noiseScale)]) , 'thisW_HomogModel','APLgains_tune_is_train',...
                    'thetaH_Ftheta_tune_is_train','thisW_Kennedy', 'thisW','theta_tune_is_train', 'thetaS_tune_is_train', 'theta_Activity_homeo_tune_is_train','InhibitionGain_tune_is_train','thisW_equalizedModel_tune_is_train','PNtrials_tune_train','theta_comp2_tune_is_train','thisW_ActivityBasedComp_tune_is_train','theta_inhibitionPlast_tune_is_train','APLgains_model6_tune_is_train','thisW_ActivityBasedComp_inhibitionPlast','thisW_ActivityBasedComp_noxjk_tuneis_train',...
                    'theta_comp2_noxjk_tuneis_train','APLgains_noxjk_tuneis_train');
                
            end
 
%         green model only
            save( strcat('robustnessTest_tunedSubgroup_trainedAllOtherGroups_BalancedEsters_and_Alcohols_tuneNo_',num2str(tune),[' _fly_wNoise',num2str(randomTrials),num2str(noiseScale),'green']) , ...
                'theta_inhibitionPlast','APLgains_model6','thisW_ActivityBasedComp_inhibitionPlast');
            
            save( strcat('robustnessTest_tuned_and_trained_onSameGroups_BalancedEsters_and_Alcohols_tuneNo_',num2str(tune),[' _fly_wNoise',num2str(randomTrials),num2str(noiseScale),'green']) , ...
                'theta_inhibitionPlast_tune_is_train','APLgains_model6_tune_is_train','thisW_ActivityBasedComp_inhibitionPlast');

            %magenta model only
            save( strcat('robustnessTest_tunedSubgroup_trainedAllOtherGroups_BalancedEsters_and_Alcohols_tuneNo_',num2str(tune),[' _fly_wNoise',num2str(randomTrials),num2str(noiseScale),'magenta']) , ...
                'thisW_Kennedy','theta_Activity_homeo','APLgains','PNtrials');
            
            save( strcat('robustnessTest_tuned_and_trained_onSameGroups_BalancedEsters_and_Alcohols_tuneNo_',num2str(tune),[' _fly_wNoise',num2str(randomTrials),num2str(noiseScale),'magenta']) , ...
                'theta_Activity_homeo_tune_is_train','thisW_Kennedy','APLgains_tune_is_train','PNtrials_tune_train');
           
            %% training and testing of the tuned models
            
            Activations=zeros(n,odors*numTrials);
            ActivationsEqualized=zeros(n,odors*numTrials);
            ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numTrials);
            Activations_comp2= zeros(n,odors*numTrials);
            Activations_theta_activity_homeo= zeros(n,odors*numTrials);
            Activations_inhibPlast=zeros(n,odors*numTrials);
            Activations_comp2_noxjk= zeros(n,odors*numTrials);
            
            
            Y=zeros(n,odors*numTrials);
            YEqualized=zeros(n,odors*numTrials);
            YHomogdummy_Ftheta=zeros(n,odors*numTrials);
            Y_comp2= zeros(n,odors*numTrials);
            Y_theta_activity_homeo=zeros(n,odors*numTrials);
            Y_inhibPlast= zeros(n,odors*numTrials);
            Y_comp2_noxjk= zeros(n,odors*numTrials);
            
            
            % KCs neural responses to odors used from training same as in tuning:
            Activations_tune_is_train=zeros(n,odors*numTrials);
            ActivationsEqualized_tune_is_train=zeros(n,odors*numTrials);
            ActivationsHomogenousdummy_Ftheta_tune_is_train=zeros(n,odors*numTrials);
            Activations_comp2_tune_is_train= zeros(n,odors*numTrials);
            Activations_theta_activity_homeo_tune_is_train= zeros(n,odors*numTrials);
            Activations_inhibPlast_tune_is_train=zeros(n,odors*numTrials);
            Activations_comp2_tune_is_train_noxjk= zeros(n,odors*numTrials);
            
            Y_tune_is_train=zeros(n,odors*numTrials);
            YEqualized_tune_is_train=zeros(n,odors*numTrials);
            YHomogdummy_Ftheta_tune_is_train=zeros(n,odors*numTrials);
            Y_comp2_tune_is_train= zeros(n,odors*numTrials);
            Y_theta_activity_homeo_tune_is_train=zeros(n,odors*numTrials);
            Y_inhibPlast_tune_is_train= zeros(n,odors*numTrials);
            Y_comp2_tune_is_train_noxjk= zeros(n,odors*numTrials);
            
            for trial = 1:(odors*numTrials)
                
                
                ActivationsHomogenousdummy_Ftheta(:,trial) = thisW_HomogModel'*PNtrials_tune_train(:,trial  );
                YHomogdummy_Ftheta(:,trial)=(( ActivationsHomogenousdummy_Ftheta(:,trial)-(APLgains(2))*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomogenousdummy_Ftheta(:,trial)-APLgains(2)*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);
                
                %
                ActivationsHomogenousdummy_Ftheta_tune_is_train(:,trial) = thisW_HomogModel'*PNtrials_tune_train(:,trial  );
                YHomogdummy_Ftheta_tune_is_train(:,trial)=(( ActivationsHomogenousdummy_Ftheta_tune_is_train(:,trial)-(APLgains_tune_is_train(2))*repmat(sum(ActivationsHomogenousdummy_Ftheta_tune_is_train(:,trial),1),n,1)-thetaH_Ftheta_tune_is_train)>0 ).*( ActivationsHomogenousdummy_Ftheta_tune_is_train(:,trial)-APLgains_tune_is_train(2)*repmat(sum(ActivationsHomogenousdummy_Ftheta_tune_is_train(:,trial),1),n,1)-thetaH_Ftheta_tune_is_train);
                
                
                
                Activations(:,trial) = thisW'*PNtrials_tune_train(:,trial );
                Y(:,trial)=(( Activations(:,trial)-(APLgains(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgains(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
                
                %
                Activations_tune_is_train(:,trial) = thisW'*PNtrials_tune_train(:,trial );
                Y_tune_is_train(:,trial)=(( Activations_tune_is_train(:,trial)-(APLgains_tune_is_train(1))*repmat(sum(Activations_tune_is_train(:,trial),1),n,1)-thetaS_tune_is_train)>0 ).*( Activations_tune_is_train(:,trial)-APLgains_tune_is_train(1)*repmat(sum(Activations_tune_is_train(:,trial),1),n,1)-thetaS_tune_is_train);
                
                
                
                ActivationsEqualized(:,trial) = thisW_equalizedModel'*PNtrials_tune_train(:,trial );
                YEqualized(:,trial)=(( ActivationsEqualized(:,trial)-(InhibitionGain)*repmat(sum(ActivationsEqualized(:,trial),1),n,1)-theta)>0 ).*( ActivationsEqualized(:,trial)-InhibitionGain*repmat(sum(ActivationsEqualized(:,trial),1),n,1)-theta);
                
                %
                ActivationsEqualized_tune_is_train(:,trial) = thisW_equalizedModel_tune_is_train'*PNtrials_tune_train(:,trial );
                YEqualized_tune_is_train(:,trial)=(( ActivationsEqualized_tune_is_train(:,trial)-(InhibitionGain_tune_is_train)*repmat(sum(ActivationsEqualized_tune_is_train(:,trial),1),n,1)-theta_tune_is_train)>0 ).*( ActivationsEqualized_tune_is_train(:,trial)-InhibitionGain_tune_is_train*repmat(sum(ActivationsEqualized_tune_is_train(:,trial),1),n,1)-theta_tune_is_train);
                
                
                Activations_comp2(:,trial) = thisW_ActivityBasedComp'*PNtrials_tune_train(:,trial );
                Y_comp2(:,trial)=(( Activations_comp2(:,trial)-(APLgains(3) )*repmat(sum(Activations_comp2(:,trial),1),n,1)-theta_comp2)>0 ).*( Activations_comp2(:,trial)-APLgains(3)*repmat(sum(Activations_comp2(:,trial),1),n,1)-theta_comp2);
                
                %
                Activations_comp2_tune_is_train(:,trial) = thisW_ActivityBasedComp_tune_is_train'*PNtrials_tune_train(:,trial );
                Y_comp2_tune_is_train(:,trial)=(( Activations_comp2_tune_is_train(:,trial)-(APLgains_tune_is_train(3) )*repmat(sum(Activations_comp2_tune_is_train(:,trial),1),n,1)-theta_comp2_tune_is_train)>0 ).*( Activations_comp2_tune_is_train(:,trial)-APLgains_tune_is_train(3)*repmat(sum(Activations_comp2_tune_is_train(:,trial),1),n,1)-theta_comp2_tune_is_train);
                
                
                Activations_theta_activity_homeo(:,trial) = thisW_Kennedy'*PNtrials_tune_train(:,trial );
                Y_theta_activity_homeo(:,trial)=(( Activations_theta_activity_homeo(:,trial)-(APLgains(4) )*repmat(sum(Activations_theta_activity_homeo(:,trial),1),n,1)-theta_Activity_homeo)>0 ).*( Activations_theta_activity_homeo(:,trial)-APLgains(4)*repmat(sum(Activations_theta_activity_homeo(:,trial),1),n,1)-theta_Activity_homeo);
                %
                Activations_theta_activity_homeo_tune_is_train(:,trial) = thisW_Kennedy'*PNtrials_tune_train(:,trial );
                Y_theta_activity_homeo_tune_is_train(:,trial)=(( Activations_theta_activity_homeo_tune_is_train(:,trial)-(APLgains_tune_is_train(4) )*repmat(sum(Activations_theta_activity_homeo_tune_is_train(:,trial),1),n,1)-theta_Activity_homeo_tune_is_train)>0 ).*( Activations_theta_activity_homeo_tune_is_train(:,trial)-APLgains_tune_is_train(4)*repmat(sum(Activations_theta_activity_homeo_tune_is_train(:,trial),1),n,1)-theta_Activity_homeo_tune_is_train);
                
                
                Activations_inhibPlast(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials_tune_train(:,trial );
                Y_inhibPlast(:,trial)=(( Activations_inhibPlast(:,trial)-(APLgains_model6 ).*repmat(sum(Activations_inhibPlast(:,trial),1),n,1)-theta_inhibitionPlast)>0 ).*( Activations_inhibPlast(:,trial)-APLgains_model6.*repmat(sum(Activations_inhibPlast(:,trial),1),n,1)-theta_inhibitionPlast);
                %
                Activations_inhibPlast_tune_is_train(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials_tune_train(:,trial );
                Y_inhibPlast_tune_is_train(:,trial)=(( Activations_inhibPlast_tune_is_train(:,trial)-(APLgains_model6_tune_is_train ).*repmat(sum(Activations_inhibPlast_tune_is_train(:,trial),1),n,1)-theta_inhibitionPlast_tune_is_train)>0 ).*( Activations_inhibPlast_tune_is_train(:,trial)-APLgains_model6_tune_is_train.*repmat(sum(Activations_inhibPlast_tune_is_train(:,trial),1),n,1)-theta_inhibitionPlast_tune_is_train);
                
                Activations_comp2_noxjk(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials_tune_train(:,trial );
                Y_comp2_noxjk(:,trial)=(( Activations_comp2_noxjk(:,trial)-(APLgains_noxjk )*repmat(sum(Activations_comp2_noxjk(:,trial),1),n,1)-theta_comp2_noxjk)>0 ).*( Activations_comp2_noxjk(:,trial)-APLgains_noxjk*repmat(sum(Activations_comp2_noxjk(:,trial),1),n,1)-theta_comp2_noxjk);
                
                %
                Activations_comp2_tune_is_train_noxjk(:,trial) = thisW_ActivityBasedComp_noxjk_tuneis_train'*PNtrials_tune_train(:,trial );
                Y_comp2_tune_is_train_noxjk(:,trial)=(( Activations_comp2_tune_is_train_noxjk(:,trial)-(APLgains_noxjk_tuneis_train )*repmat(sum(Activations_comp2_tune_is_train_noxjk(:,trial),1),n,1)-theta_comp2_noxjk_tuneis_train)>0 ).*( Activations_comp2_tune_is_train_noxjk(:,trial)-APLgains_noxjk_tuneis_train*repmat(sum(Activations_comp2_tune_is_train_noxjk(:,trial),1),n,1)-theta_comp2_noxjk_tuneis_train);
                
                
            end
            
            %% training of the KC-MBONs output weights
            
            for l_r=1:lrs
                
                WopAllOdours=1*rand(n,2);
                WopAllOdoursEqualized= WopAllOdours;
                WopAllOdoursHomog_Ftheta=WopAllOdours;
                WopAllOdoursInhPlast=WopAllOdours;
                WopAllOdoursEqualizedComp2=WopAllOdours;
                WopFromPNs= 1*rand(m,2);
                WopAllOdoursThetaActivityHomeo=WopAllOdours;
                WopAllOdoursEqualizedComp2_noxjk=WopAllOdours;
                
                % the weights vectors initial state in the tune_is_train
                % models should be equal to those with tuning set different
                % than training set
                
                WopAllOdours_tune_is_train=WopAllOdours;
                WopAllOdoursEqualized_tune_is_train= WopAllOdours;
                WopAllOdoursHomog_Ftheta_tune_is_train=WopAllOdours;
                WopAllOdoursInhPlast_tune_is_train=WopAllOdours;
                WopAllOdoursEqualizedComp2_tune_is_train=WopAllOdours;
                WopFromPNs_tune_is_train= 1*rand(m,2);
                WopAllOdoursThetaActivityHomeo_tune_is_train=WopAllOdours;
                WopAllOdoursEqualizedComp2_tune_is_train_noxjk=WopAllOdours;
                
                
                alpha=0.000001* (10^((l_r)));
                
                c=1;
                
                YHomog_Fthetatemp=reshape(YHomogdummy_Ftheta,n,odors,numTrials);
                
                Ytemp= reshape(Y,n,odors,numTrials);
                
                YEqualizedtemp=reshape(YEqualized,n,odors,numTrials);
                
                Y_comp2temp= reshape(Y_comp2,n,odors,numTrials);
                
                Y_theta_activity_homeotemp= reshape(Y_theta_activity_homeo,n,odors,numTrials);
                
                Y_inhibPlasttemp= reshape(Y_inhibPlast,n,odors,numTrials);
                Y_comp2temp_noxjk= reshape(Y_comp2_noxjk,n,odors,numTrials);
                
                
                
                % rescale the KC responses to be from 0-1
                Ytemp=rescale(Ytemp);
                YEqualizedtemp=rescale(YEqualizedtemp);
                YHomog_Fthetatemp=rescale(YHomog_Fthetatemp);
                Y_comp2temp=rescale(Y_comp2temp);
                Y_comp2temp_noxjk=rescale(Y_comp2temp_noxjk);
                Y_theta_activity_homeotemp= rescale(Y_theta_activity_homeotemp);
                Y_inhibPlasttemp=rescale(Y_inhibPlasttemp);
                
                YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
                
                Ytr=Ytemp(:,:,1:numtrainingSamples);
                
                YEqualizedtr=YEqualizedtemp(:,:,1:numtrainingSamples);
                
                Y_comp2tr= Y_comp2temp(:,:,1:numtrainingSamples);
                
                Y_theta_activity_homeotr= Y_theta_activity_homeotemp(:,:,1:numtrainingSamples);
                
                Y_inhibPlasttr= Y_inhibPlasttemp(:,:,1:numtrainingSamples);
                Y_comp2tr_noxjk= Y_comp2temp_noxjk(:,:,1:numtrainingSamples);
                
                
                % rescaling and reshaping the responses of the
                % tune_is_train data
                
                YHomog_Fthetatemp_tune_is_train=reshape(YHomogdummy_Ftheta_tune_is_train,n,odors,numTrials);
                
                Ytemp_tune_is_train= reshape(Y_tune_is_train,n,odors,numTrials);
                
                YEqualizedtemp_tune_is_train=reshape(YEqualized_tune_is_train,n,odors,numTrials);
                
                Y_comp2temp_tune_is_train= reshape(Y_comp2_tune_is_train,n,odors,numTrials);
                
                Y_theta_activity_homeotemp_tune_is_train= reshape(Y_theta_activity_homeo_tune_is_train,n,odors,numTrials);
                
                Y_inhibPlasttemp_tune_is_train= reshape(Y_inhibPlast_tune_is_train,n,odors,numTrials);
                Y_comp2temp_tune_is_train_noxjk= reshape(Y_comp2_tune_is_train_noxjk,n,odors,numTrials);
                
                
                
                % rescale the KC responses to be from 0-1
                Ytemp_tune_is_train=rescale(Ytemp_tune_is_train);
                
                YEqualizedtemp_tune_is_train=rescale(YEqualizedtemp_tune_is_train);
                
                YHomog_Fthetatemp_tune_is_train=rescale(YHomog_Fthetatemp_tune_is_train);
                
                Y_comp2temp_tune_is_train=rescale(Y_comp2temp_tune_is_train);
                
                Y_theta_activity_homeotemp_tune_is_train= rescale(Y_theta_activity_homeotemp_tune_is_train);
                
                Y_inhibPlasttemp_tune_is_train=rescale(Y_inhibPlasttemp_tune_is_train);
                
                Y_comp2temp_tune_is_train_noxjk=rescale(Y_comp2temp_tune_is_train_noxjk);
                
                
                YHomog_Fthetatr_tune_is_train=YHomog_Fthetatemp_tune_is_train(:,:,1:numtrainingSamples);
                
                Ytr_tune_is_train=Ytemp_tune_is_train(:,:,1:numtrainingSamples);
                
                YEqualizedtr_tune_is_train=YEqualizedtemp_tune_is_train(:,:,1:numtrainingSamples);
                
                Y_comp2tr_tune_is_train= Y_comp2temp_tune_is_train(:,:,1:numtrainingSamples);
                
                Y_theta_activity_homeotr_tune_is_train= Y_theta_activity_homeotemp_tune_is_train(:,:,1:numtrainingSamples);
                
                Y_inhibPlasttr_tune_is_train= Y_inhibPlasttemp_tune_is_train(:,:,1:numtrainingSamples);
                
                Y_comp2tr_tune_is_train_noxjk= Y_comp2temp_tune_is_train_noxjk(:,:,1:numtrainingSamples);
                
                %% learning from a preceptron in the output layer; by Synaptic depression of the
                % input weights to the MBON of the opposing valence.
                
                for odour=1:odors
                    
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
                        
                        delta = exp( -(alpha/mean(YEqualizedtr(:)))* sum((YEqualizedtr(:,odour,:)),3) );
                        
                        
                        WopAllOdoursEqualized(:,2)= WopAllOdoursEqualized(:,2).*delta;
                        
                    else
                        
                        delta =  exp(-(alpha /mean(YEqualizedtr(:)))* sum(YEqualizedtr(:,odour,:),3) );
                        WopAllOdoursEqualized(:,1)= WopAllOdoursEqualized(:,1) .*delta;
                        
                    end
                    
                    
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_comp2tr(:)))* sum(Y_comp2tr(:,odour,:),3) );
                        
                        WopAllOdoursEqualizedComp2(:,2)= WopAllOdoursEqualizedComp2(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_comp2tr(:)))* sum(Y_comp2tr(:,odour,:),3) );
                        WopAllOdoursEqualizedComp2(:,1)= WopAllOdoursEqualizedComp2(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_theta_activity_homeotr(:)))* sum(Y_theta_activity_homeotr(:,odour,:),3) );
                        
                        WopAllOdoursThetaActivityHomeo(:,2)= WopAllOdoursThetaActivityHomeo(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_theta_activity_homeotr(:)))* sum(Y_theta_activity_homeotr(:,odour,:),3) );
                        WopAllOdoursThetaActivityHomeo(:,1)= WopAllOdoursThetaActivityHomeo(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_inhibPlasttr(:)))* sum(Y_inhibPlasttr(:,odour,:),3) );
                        
                        WopAllOdoursInhPlast(:,2)= WopAllOdoursInhPlast(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_inhibPlasttr(:)))* sum(Y_inhibPlasttr(:,odour,:),3) );
                        WopAllOdoursInhPlast(:,1)= WopAllOdoursInhPlast(:,1) .*delta;
                        
                    end
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_comp2tr_noxjk(:)))* sum(Y_comp2tr_noxjk(:,odour,:),3) );
                        
                        WopAllOdoursEqualizedComp2_noxjk(:,2)= WopAllOdoursEqualizedComp2_noxjk(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_comp2tr_noxjk(:)))* sum(Y_comp2tr_noxjk(:,odour,:),3) );
                        WopAllOdoursEqualizedComp2_noxjk(:,1)= WopAllOdoursEqualizedComp2_noxjk(:,1) .*delta;
                        
                    end
                    
                    
                    %% do training for the output perceptrons from the tune_is_train responses
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        delta =  exp( -(alpha/mean(YHomog_Fthetatr_tune_is_train(:)))* sum(YHomog_Fthetatr_tune_is_train(:,odour,:),3) );
                        
                        WopAllOdoursHomog_Ftheta_tune_is_train(:,2)= WopAllOdoursHomog_Ftheta_tune_is_train(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(YHomog_Fthetatr_tune_is_train(:)))* sum(YHomog_Fthetatr_tune_is_train(:,odour,:),3) );
                        WopAllOdoursHomog_Ftheta_tune_is_train(:,1)= WopAllOdoursHomog_Ftheta_tune_is_train(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        delta =  exp( -(alpha/mean(Ytr_tune_is_train(:)))* sum(Ytr_tune_is_train(:,odour,:),3) );
                        
                        WopAllOdours_tune_is_train(:,2)= WopAllOdours_tune_is_train(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Ytr_tune_is_train(:)))* sum(Ytr_tune_is_train(:,odour,:),3) );
                        WopAllOdours_tune_is_train(:,1)= WopAllOdours_tune_is_train(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta = exp( -(alpha/mean(YEqualizedtr_tune_is_train(:)))* sum((YEqualizedtr_tune_is_train(:,odour,:)),3) );
                        
                        
                        WopAllOdoursEqualized_tune_is_train(:,2)= WopAllOdoursEqualized_tune_is_train(:,2).*delta;
                        
                    else
                        
                        delta =  exp(-(alpha /mean(YEqualizedtr_tune_is_train(:)))* sum(YEqualizedtr_tune_is_train(:,odour,:),3) );
                        WopAllOdoursEqualized_tune_is_train(:,1)= WopAllOdoursEqualized_tune_is_train(:,1) .*delta;
                        
                    end
                    
                    
                    
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_comp2tr_tune_is_train(:)))* sum(Y_comp2tr_tune_is_train(:,odour,:),3) );
                        
                        WopAllOdoursEqualizedComp2_tune_is_train(:,2)= WopAllOdoursEqualizedComp2_tune_is_train(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_comp2tr_tune_is_train(:)))* sum(Y_comp2tr_tune_is_train(:,odour,:),3) );
                        WopAllOdoursEqualizedComp2_tune_is_train(:,1)= WopAllOdoursEqualizedComp2_tune_is_train(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_theta_activity_homeotr_tune_is_train(:)))* sum(Y_theta_activity_homeotr_tune_is_train(:,odour,:),3) );
                        
                        WopAllOdoursThetaActivityHomeo_tune_is_train(:,2)= WopAllOdoursThetaActivityHomeo_tune_is_train(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_theta_activity_homeotr_tune_is_train(:)))* sum(Y_theta_activity_homeotr_tune_is_train(:,odour,:),3) );
                        
                        WopAllOdoursThetaActivityHomeo_tune_is_train(:,1)= WopAllOdoursThetaActivityHomeo_tune_is_train(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_inhibPlasttr_tune_is_train(:)))* sum(Y_inhibPlasttr_tune_is_train(:,odour,:),3) );
                        
                        WopAllOdoursInhPlast_tune_is_train(:,2)= WopAllOdoursInhPlast_tune_is_train(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_inhibPlasttr_tune_is_train(:)))* sum(Y_inhibPlasttr_tune_is_train(:,odour,:),3) );
                        WopAllOdoursInhPlast_tune_is_train(:,1)= WopAllOdoursInhPlast_tune_is_train(:,1) .*delta;
                        
                    end
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_comp2tr_tune_is_train_noxjk(:)))* sum(Y_comp2tr_tune_is_train_noxjk(:,odour,:),3) );
                        
                        WopAllOdoursEqualizedComp2_tune_is_train_noxjk(:,2)= WopAllOdoursEqualizedComp2_tune_is_train_noxjk(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_comp2tr_tune_is_train_noxjk(:)))* sum(Y_comp2tr_tune_is_train_noxjk(:,odour,:),3) );
                        WopAllOdoursEqualizedComp2_tune_is_train_noxjk(:,1)= WopAllOdoursEqualizedComp2_tune_is_train_noxjk(:,1) .*delta;
                        
                    end
                    
                    
                    
                end
                
                %% perfromance as a function of the strictness of the decision making
                %% this strictness is dictated by C in the soft-max function.
                %% so given the same fly, same task, and after learning measure the performance as f(c)
                
                for c=1:Crange
                    
                    
                    C=C_SoftMax*(10^c);
                    
                    [acc,accEq,accEq2, accKenn,accInhPlast,acc_comp2_noxjk]=Robustness_testing(C,WopAllOdours,WopAllOdoursEqualized,WopAllOdoursEqualizedComp2, WopAllOdoursThetaActivityHomeo, WopAllOdoursInhPlast, WopAllOdoursEqualizedComp2_noxjk,...
                        classAction1,numTrials,numtrainingSamples,Ytemp,YEqualizedtemp,Y_comp2temp,Y_theta_activity_homeotemp,Y_inhibPlasttemp,Y_comp2temp_noxjk);
                    
                    test_p_raEq(randomTrials,tune,l_r,c)=accEq;
                    test_p_ra(randomTrials,tune,l_r,c)=acc;
                    test_p_raEq2(randomTrials,tune,l_r,c)=accEq2;
                    test_p_raThetaActivity_homeo(randomTrials,tune,l_r,c)= accKenn;
                    test_p_raInhPlast(randomTrials,tune,l_r,c)= accInhPlast;
                    test_p_raEq2_noxjk(randomTrials,tune,l_r,c)=acc_comp2_noxjk;
                    
                    
                    % testing of the tuned models in a familiar environment.
                    [acc_tune_is_train,accEq_tune_is_train,accEq2_tune_is_train, accKenn_tune_is_train,accInhPlast_tune_is_train,accEq2_tune_is_train_noxjk]=Robustness_testing(C,WopAllOdours_tune_is_train,WopAllOdoursEqualized_tune_is_train,WopAllOdoursEqualizedComp2_tune_is_train, WopAllOdoursThetaActivityHomeo_tune_is_train, WopAllOdoursInhPlast_tune_is_train,WopAllOdoursEqualizedComp2_tune_is_train_noxjk,...
                        classAction1,numTrials,numtrainingSamples,Ytemp_tune_is_train,YEqualizedtemp_tune_is_train,Y_comp2temp_tune_is_train,Y_theta_activity_homeotemp_tune_is_train,Y_inhibPlasttemp_tune_is_train,Y_comp2temp_tune_is_train_noxjk);
                    
                    test_p_raEq_tune_is_train(randomTrials,tune,l_r,c)=accEq_tune_is_train;
                    test_p_ra_tune_is_train(randomTrials,tune,l_r,c)=acc_tune_is_train;
                    test_p_raEq2_tune_is_train(randomTrials,tune,l_r,c)=accEq2_tune_is_train;
                    test_p_raThetaActivity_homeo_tune_is_train(randomTrials,tune,l_r,c)= accKenn_tune_is_train;
                    test_p_raInhPlast_tune_is_train(randomTrials,tune,l_r,c)= accInhPlast_tune_is_train;
                    test_p_raEq2_noxjk_tune_is_train(randomTrials,tune,l_r,c)=accEq2_tune_is_train_noxjk;
                    
                    
                    
                    [accH1]=HomogenousModel_KernelTesting(C,WopAllOdoursHomog_Ftheta,PNtrials_tune_train, thetaH_Ftheta,InhibitionGain, APLgains,classAction1,numTrials,numtrainingSamples,YHomog_Fthetatemp);
                    
                    test_p_raH_FixedTheta(randomTrials,tune,l_r,c)=accH1;
                    
                    
                    [accH1_tune_is_train]=HomogenousModel_KernelTesting(C,WopAllOdoursHomog_Ftheta_tune_is_train,PNtrials_tune_train, thetaH_Ftheta_tune_is_train,InhibitionGain_tune_is_train, APLgains_tune_is_train,classAction1,numTrials,numtrainingSamples,YHomog_Fthetatemp_tune_is_train);
                    
                    test_p_raH_FixedTheta_tune_is_train(randomTrials,tune,l_r,c)=accH1_tune_is_train;
                    
                end
                
                
                
            end
            
            
        end
        
     disp (strcat('tune',num2str(tune)))   
    end
end

%% saving the model accuracies in the novel and familiar environments
%% for all the models, in each tune for a certain chemical group
save('test_p_raEq.mat','test_p_raEq');
save('test_p_raEq_tune_is_train.mat','test_p_raEq_tune_is_train')
save('test_p_raRand.mat','test_p_ra')
save('test_p_raRand_tune_is_train.mat','test_p_ra_tune_is_train')
save('test_p_ra_H_W.mat','test_p_raEq2')
save('test_p_ra_H_W_tune_is_train.mat','test_p_raEq2_tune_is_train')
save('test_p_ra_H_W_noxjk.mat','test_p_raEq2_noxjk')
save('test_p_ra_H_W_noxjk_tune_is_train.mat','test_p_raEq2_noxjk_tune_is_train')
save('test_p_ra_H_Theta.mat','test_p_raThetaActivity_homeo')
save('test_p_ra_H_Theta_tune_is_train.mat','test_p_raThetaActivity_homeo_tune_is_train')
save('test_p_ra_H_InhPlast.mat','test_p_raInhPlast')
save('test_p_ra_H_InhPlast_tune_is_train.mat','test_p_raInhPlast_tune_is_train')
save('test_p_raHomog_FixedTheta.mat','test_p_raH_FixedTheta')
save('test_p_raHomog_FixedTheta_tune_is_train.mat','test_p_raH_FixedTheta_tune_is_train')

%% END  OF CODE------------------------------------------------------------
delete(gcp('nocreate'));