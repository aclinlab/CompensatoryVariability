%% this code replicates the results for the main compensatory models in Fig 5
% it calculates the performances scroes, Lifetime sparsity and their
% dimensionality.
% it also does the tuning and calculates the performance
%for the alternative rules in the supplementary materials:
% the alternative models for tuning input PN-KC synaptic weights: FigS3
% the alternative model of tuning theta for equalizing KCs firing probabilities
% (Kennedy's Inspired model)

%% Start of the code

clc
clear all;

CurrDir= pwd;

load('hallem_olsen.mat');
% load('PW_given_N.mat');
% load('PW_given_theta_and_n.mat');
% load('W_PN_KC.mat');
% load ('P_n.mat');
% ll=20;
ll=1;
numtrainingSamples=15; %% artificial plus linearly dependent responses
modss=1;
mulOd=110;
lrs=7;
Crange=1;
odors=mulOd;
numTrials = 30;
PN = hallem_olsen(1:110,:)';
PNs=zeros(24,odors);
theta_Vs_deltaActivations=[ 0, 0];
k=0;
C_SoftMax=1;
NScales=1;

%% create artificial odors, n odors, drawn from the PNs responses distribution
% in Hallem-Carlson data:
for Pn=1:24
    [prob,bins]=hist(PN(Pn,:),100);
    prob=prob/sum(prob);

    PNs(Pn,k+1:k+odors)=randsample(bins,odors,'true',prob);

end

%% start a local parallel pool of 3 threads to speed up simulation time
p=parpool(3);
p.IdleTimeout = 1000;
parfevalOnAll(@maxNumCompThreads,0,3)


% for mods=1:modss
mods=1
    
    odors= mulOd*(mods/modss);
    x=PN;
    
    for randomTrials=1:ll
        
        classAction1=randsample([1:odors],round(odors/2),'false');
        
        % theta_threshold= 0.0+(10);
        %% HO data IG=0.6
        
        InhibitionGain= 0.0+ (0.0);
        
        APLgainP= zeros(1,2);
        % n = MBmodel.nKCs; % number of neurons in the hidden layer
        % m = MBmodel.nPNs;  %number of dimensions in the input data
        n=2000;
        m=24;
        
        
        % randomly assign thresholds for KCs
        % ThetaMu=13;
        % theta=abs(normrnd(ThetaMu,ThetaMu*(5.6/21.5),[n 1])); %% avoid negative values of theta
        % theta(theta>70)=70;
        % theta(theta<0.01)=0.01;
        
        %initialize the weights matrix between KC and PN
        thisW = zeros(m, n);
        % thisW_equalizedModel=zeros(m,n);
        % thisW_HomogModel=zeros(m,n);
        % thisW_Kennedy= zeros(m,n);
        % initial_thisW_ActivityBasedComp= zeros(m,n);
        
        % default model, a.k.a. random model, all variable, without tuning
        MBmodel_varN_varW_varTh = MBmodelBuilder(2000,24,true,true,true, 'theta_mean',10, 'theta_std',10*5.6/21.5 );
        % Homog model: fix Nclaws, fix weights, 
        % MBmodel_fixN_fixW_varTh = MBmodelBuilder(2000,24,false,false,false, 'theta_mean',5+randn(1), 'theta_std',5*5.6/21.5 );
        % Kennedy model: var Nclaws, var weights, 
        % equalized model: var Nclaws, var weights(N,theta), 
        % activityComp: var Nclaws, var weight=5+rand
        
        mComp=0.1;
        constraints=0;
        
        
        for noiseScale=1:NScales
            
            noise=((noiseScale^2)/10);
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
            
            PNtrials(PNtrials<0)=0;
            PNtrials=rescale(PNtrials,0,5);
            
            tic
            C_theta=1;
            
            spmd
                
                ActivationsDummy=zeros(n,odors*numtrainingSamples);
                ActivationsEqualizeddummy=zeros(n,odors*numtrainingSamples);
                A=zeros(n,odors*numtrainingSamples);
                Y=zeros(n,odors*numtrainingSamples);
                YEqualizeddummy=zeros(n,odors*numtrainingSamples);
                codingLevelEqualizedDummy=[];
                
                if labindex==1
                    % tuning the activity independent compensatory model
                    % P(W|n, theta)
                    InhAbs_CLVec(1)=0;
                    InhAbsTarget=0.20;

                    while(~constraints)
                        
                        % with inhibition gain absent, scaling thetas distribution to acheive
                        % the sparsity constraint= 20% when APL feedback is inhibited.
                        
                        eta=2;
                        for trial = 1:(odors*numtrainingSamples)
                            
                            A(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_theta.*theta))>0 ).*( A(:,trial)-(C_theta.*theta));
                            codingLevelEqualizedDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        InhAbs_mComp=mean(codingLevelEqualizedDummy);
                        depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                        depsi1_dy(isnan(depsi1_dy))=0;
                        depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(theta,1,odors*numtrainingSamples));
                        Grad= ((InhAbs_mComp)-InhAbsTarget)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                        C_theta=C_theta -(eta.*Grad);
                        
                        if (C_theta<0)
                            error('the scale factor in cyan model is -ve')
                        end
                        
                        %% resample the weights
                        %%  ---> WHY??
                        thisW_equalizedModel=zeros(m,n);
                        
                        for k=1:n
                            PNpartners = PnToKc{k};
                            for j=1:length(PNpartners)
                                whichPN = PNpartners(j);
                                % pick random weight from a log normal distribution that
                                % roughtly fits the Turner distribution
                                % sample the weights from the fitted weights in the other script (modelling KC_PnWeights.m)
                                
                                ThetaInd= round(((theta(k)-0.01)/0.1)+1);
                                % capping theta at 0.1..
                                ThetaInd(ThetaInd==0)=1;
                                this_KCWeights= PW_given_theta_and_n(length(PnToKc{k})-1,ThetaInd,:);
                                thisWeight_equalizedModel= randsample(W,1,'true', this_KCWeights);
                                thisW_equalizedModel(whichPN,k)= thisW_equalizedModel(whichPN,k)+thisWeight_equalizedModel;
                            end
                        end
                        
                        % replicating the sparsity level of the KCs in real MB network
                        % CL=10%
                        eta_2=0.000001;
                        for trial = 1:(odors*numtrainingSamples)
                            ActivationsEqualizeddummy(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            YEqualizeddummy(:,trial)=(( ActivationsEqualizeddummy(:,trial)-(InhibitionGain)*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta.*theta))>0 ).*( ActivationsEqualizeddummy(:,trial)-InhibitionGain*repmat(sum(ActivationsEqualizeddummy(:,trial),1),n,1)-(C_theta.*theta));
                            codingLevelEqualizedDummy(trial)=  (sum(YEqualizeddummy(:,trial)>0,1)/n);
                        end
                        mComp=mean(codingLevelEqualizedDummy);
                        
                        dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                        dsig_dy(isnan(dsig_dy))=0;
                        dAct_dalpha= sum(ActivationsEqualizeddummy,1);
                        dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                        Grad_alpha= ((mComp)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                        InhibitionGain= InhibitionGain- eta_2*(Grad_alpha);
                        
                        if (InhibitionGain<0)
                            wh=1;
                        end
                        
                        % check if the constraints are satisfied
                        
                        for trial = 1:(odors*numtrainingSamples)
                            
                            ActivationsDummy(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_theta.*theta))>0 ).*( ActivationsDummy(:,trial)-(C_theta.*theta));
                            codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                        end
                        
                        InhAbs_CL=mean(codingLevelDummy);
                        
                        for trial = 1:(odors*numtrainingSamples)
                            
                            Activations(:,trial) = thisW_equalizedModel'*PNtrials(:,trial);
                            Y(:,trial)=(( Activations(:,trial)-(InhibitionGain)*repmat(sum(Activations(:,trial),1),n,1)-(C_theta.*theta))>0 ).*( Activations(:,trial)-InhibitionGain*repmat(sum(Activations(:,trial),1),n,1)-(C_theta.*theta));
                            codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);
                            
                        end
                        CL_=mean(codingLevelDummy);
                        constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.2) &( abs(CL_-0.10)<0.01 );
                        InhAbs_CLVec(end+1)=InhAbs_CL;
                        
                        % debugging
                        theta_Vs_deltaActivations(end+1,:)=[ mean(theta), mean(ActivationsDummy(find(ActivationsDummy)))];
                        
                        
                    end % end of optimisation
                    theta=(C_theta.*theta);
                    CLevelP(1)=CL_;
                    INHAbs_CLP(1)=InhAbs_CL;
                    
                    
                elseif labindex==2
                    MBmodel_varN_varW_varTh = optimiseMBparams(MBmodel_varN_varW_varTh, reshape(PNtrials,size(PNtrials,1),[]) );
                elseif labindex==3
                    % replicating the KCs sparsity levels in the homogenous model
                    % no tuning of any activity dependent plasticity.
                    MBmodel_fixN_fixW_fixTh = optimiseMBparams(MBmodel_fixN_fixW_fixTh, reshape(PNtrials,size(PNtrials,1),[]) );
                end
            end
            Clevels=[CLevelP{1}, MBmodel_varN_varW_varTh.CL_incInh , MBmodel_fixN_fixW_fixTh.CL_incInh];
            INHAbs_CL=[INHAbs_CLP{1}, MBmodel_varN_varW_varTh.CL_disInh , MBmodel_fixN_fixW_fixTh.CL_disInh];
            
            % end of the parallel pool, aggregating the tuned variables of the KCs
            % from the 3 jobs/threads into a matrix
            % for i=1:2
            %     temp=APLgainP{i+1};
            %     APLgains(i)= temp(i);
            % end
            % replaced by
            APLgains(2) = MBmodel_varN_varW_varTh.alpha;
            APLgains(3) = MBmodel_fixN_fixW_varTh.alpha;
            
            InhibitionGain=InhibitionGain{1};
            theta=theta{1};
            thetaH_Ftheta=thetaH_Ftheta{3};
            thetaS=thetaS{2};
            
            thisW_equalizedModel=thisW_equalizedModel{1};
            
            disp('time to run black,red,cyan models:');
            toc
            
            tic
            %% tuning the input PN-KC weights using the simplified learning rule
            % our main model for tuning the PN-KC weights presented in Fig5 (the Blue
            % % model)
            thisW_ActivityBasedComp_noxjk=initial_thisW_ActivityBasedComp;
            Activations =zeros(n,odors*numtrainingSamples);
            ActivationsDummy= zeros(n, odors*numtrainingSamples);
            Conn=zeros(m,n);
            Conn(find(thisW_ActivityBasedComp_noxjk))=1;
            mask= zeros (m,n);
            mask(find(thisW_ActivityBasedComp_noxjk))=1;
            C_=1;
            APLgains_noxjk=0;
            T=5;
            A0=(0.51)*ones(n,1);
            epsilon= A0(1)*0.06;
            theta_comp2_noxjk= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_comp2_noxjk(theta_comp2_noxjk>70)=70;
            theta_comp2_noxjk(theta_comp2_noxjk<0.01)=0.01;
            Inhabs_CLV=[];
            CLV=[];
            conditions=0; %%
            A=zeros(n,odors*numtrainingSamples);
            Y_d=zeros(n,odors*numtrainingSamples);
            codingLevelDummy=[];
            t=1;
            while(~conditions)
                
                % with inhibition gain absent, scaling thetas distribution to acheive
                % the sparsity constraint= 20% when APL feedback is inhibited.
                % adjust C_theta
                for trial = 1:(odors*numtrainingSamples)
                    A(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_.*theta_comp2_noxjk) )>0 ).*( A(:,trial)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end                
                InhAbs_CL=mean(codingLevelDummy);
                depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_comp2_noxjk,1,odors*numtrainingSamples));
                %eta=10.*mean(Y_d,2);
                eta=10;
                Grad= ((InhAbs_CL)-0.20)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
                C_=C_ - (eta*Grad);
                if (C_<0)
                    error('the scale factor in the blue model is -ve!!')
                end
                
                % replicating the sparsity level of the KCs in real MB network
                % CL=10%
                % adjust APLgain
                eta_2=0.00000005;
                for trial = 1:(odors*numtrainingSamples)
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
                Grad_alpha= ((CL_)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains_noxjk= APLgains_noxjk- eta_2*(Grad_alpha);
                
                %now, third optimization step finding the optimum weights
                %to equalize the average activity levels among KCs, all
                %synapses are scaled by the same factor, homogenously.
                
                for trial = 1:(odors*numtrainingSamples)
                    Activations(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial );
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk))>0 ).*( Activations(:,trial)-APLgains_noxjk*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk));
                end
                avgAKcs=mean(Y_,2);
                errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);
                thisW_ActivityBasedComp_noxjk= thisW_ActivityBasedComp_noxjk-(0.05).*((1.*(mask.*errorInActivity)));
                thisW_ActivityBasedComp_noxjk(find(thisW_ActivityBasedComp_noxjk<0))=0;
                % check if the constraints are statsified:
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    ActivationsDummy(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_comp2_noxjk))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    Activations(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk))>0 ).*( Activations(:,trial)-APLgains_noxjk*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                if mod(t,20)==0
                    disp(length(find(abs(avgAKcs-A0)>epsilon)))
                    disp(CL_)
                end
                t=t+1;

                conditions= all(abs(avgAKcs-A0)<epsilon) &( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01 );
            end
            
            theta_comp2_noxjk=(C_.*theta_comp2_noxjk);
            disp('time to run blue model');
            toc
            
            
            %% tuning input excitatory weights, H_indiv in FigS3, the dark blue model
            % each PN-KC is tuned individually using the extra <x_jk>k factor
            %                 tic
            
            %                 thisW_ActivityBasedComp= initial_thisW_ActivityBasedComp;
            %                 Activations =zeros(n,odors*numtrainingSamples);
            %                 ActivationsDummy= zeros(n, odors*numtrainingSamples);
            %                 Conn=zeros(m,n);
            %                 Conn(find(thisW_ActivityBasedComp))=1;
            %                 C_=1;
            %                 APLgains(3)=0;
            %                 T=5;
            %                 A0=(0.51)*ones(n,1);
            %                 epsilon= A0(1)*0.06;
            %                 theta_comp2= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            %                 theta_comp2(theta_comp2>70)=70;
            %                 theta_comp2(theta_comp2<0.01)=0.01;
            %                 Inhabs_CLV=[];
            %                 CLV=[];
            %                 conditions=0; %%
            %                 A=zeros(n,odors*numtrainingSamples);
            %                 Y_d=zeros(n,odors*numtrainingSamples);
            %                 codingLevelDummy=[];
            %
            %                 while(~conditions)
            %
            %
            %                 % with inhibition gain absent, scaling thetas distribution to acheive
            %                 % the sparsity constraint= 20% when APL feedback is inhibited.
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     A(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
            %                     Y_d(:,trial)=(( A(:,trial)-(C_.*theta_comp2) )>0 ).*( A(:,trial)-(C_.*theta_comp2));
            %                     codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
            %                 end
            %
            %                 InhAbs_CL=mean(codingLevelDummy);
            %                 depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
            %                 depsi1_dy(isnan(depsi1_dy))=0;
            %                 depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_comp2,1,odors*numtrainingSamples));
            %                 eta=10;
            %                 Grad= ((InhAbs_CL)-0.20)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
            %                 C_=C_ - (eta*Grad);
            %
            %                 if (C_<0)
            %                     error('the scale factor in the blue model is -ve!!')
            %                 end
            %
            %                 % replicating the sparsity level of the KCs in real MB network
            %                 % CL=10%
            %                 eta_2=0.00000005;
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     Activations(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
            %                     Y_(:,trial)=(( Activations(:,trial)-(APLgains(3))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2))>0 ).*( Activations(:,trial)-APLgains(3)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2));
            %                     codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
            %
            %                 end
            %                 CL_=mean(codingLevelDummy);
            %                 dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
            %                 dsig_dy(isnan(dsig_dy))=0;
            %                 dAct_dalpha= sum(Activations,1);
            %
            %                 if ( any(isinf(dAct_dalpha)) )
            %                 stopp=1;
            %                 end
            %                 dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
            %                 Grad_alpha= ((CL_)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
            %                 APLgains(3)= APLgains(3)- eta_2*(Grad_alpha);
            %
            %                 %now, third optimization step finding the optimum weights
            %                 %to equalize the average activity levels among KCs,
            %                 %including the <x_jk>_k factor
            %
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     Activations(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial );
            %                     Y_(:,trial)=(( Activations(:,trial)-(APLgains(3))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2))>0 ).*( Activations(:,trial)-APLgains(3)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2));
            %
            %                 end
            %                 avgAKcs=mean(Y_,2);
            %                 errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);
            %                 Conn2=repmat(Conn,1,1,odors*numtrainingSamples);
            %                 Xik=reshape(PNtrials(:,1:odors*numtrainingSamples),m,1,odors*numtrainingSamples);
            %                 Xik_=repmat(Xik,1,n,1);
            %                 filt_=  ( ((1-(InhibitionGain/n)) .* (Xik_.* Conn2)) ) ;
            %                 filt_Xjm= mean(filt_,3);
            %                 mask=filt_Xjm;
            %                 thisW_ActivityBasedComp= thisW_ActivityBasedComp-(0.1).*((1.*(mask.*errorInActivity)));
            %
            %                 if (~isempty(find(isinf(thisW_ActivityBasedComp) )))
            %                 g=1;
            %                 end
            %                 % make sure all PNs-KCs weights are +ve or zero.
            %                 thisW_ActivityBasedComp(find(thisW_ActivityBasedComp<0))=0;
            %
            %                 % checking the constraints if satisfied:
            %
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     ActivationsDummy(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
            %                     Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_comp2))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_comp2));
            %                     codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
            %                 end
            %
            %                 InhAbs_CL=mean(codingLevelDummy);
            %
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     Activations(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial);
            %                     Y_(:,trial)=(( Activations(:,trial)-(APLgains(3))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2))>0 ).*( Activations(:,trial)-APLgains(3)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2));
            %                     codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
            %
            %                 end
            %                 CL_=mean(codingLevelDummy);
            %
            %                 conditions= all(abs(avgAKcs-A0)<epsilon) &( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01 );
            %
            %                 Inhabs_CLV(end+1)=InhAbs_CL;
            %                 CLV(end+1)=CL_;
            %
            %                 end
            %
            %                 theta_comp2=(C_.*theta_comp2);
            %                 toc
            %% tuning the input PN-KC weights using the learning rule derived
            %  from the error function: with the H(Y) term; H_active in FigS3
            
            %                 thisW_ActivityBasedComp_wHy= initial_thisW_ActivityBasedComp;
            %                 Activations =zeros(n,odors*numtrainingSamples);
            %                 ActivationsDummy= zeros(n, odors*numtrainingSamples);
            %                 Conn=zeros(m,n);
            %                 Conn(find(thisW_ActivityBasedComp_wHy))=1;
            %                 mask= zeros (m,n);
            %                 mask(find(thisW_ActivityBasedComp_wHy))=1;
            %                 C_=1;
            %                 APLgains_wHy=0;
            %                 T=5;
            %                 A0=(0.51)*ones(n,1);
            %                 epsilon= A0(1)*0.06;
            %                 theta_comp2_wHy= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            %                 theta_comp2_wHy(theta_comp2_wHy>70)=70;
            %                 theta_comp2_wHy(theta_comp2_wHy<0.01)=0.01;
            %                 Inhabs_CLV=[];
            %                 CLV=[];
            %                 conditions=0; %%
            %                 A=zeros(n,odors*numtrainingSamples);
            %                 Y_d=zeros(n,odors*numtrainingSamples);
            %                 Y_=zeros(n,odors*numtrainingSamples);
            %                 codingLevelDummy=[];
            %
            %                 while(~conditions)
            %
            %
            %                 % with inhibition gain absent, scaling thetas distribution to acheive
            %                 % the sparsity constraint= 20% when APL feedback is inhibited.
            %
            %                        for trial = 1:(odors*numtrainingSamples)
            %
            %                         A(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial);
            %                         Y_d(:,trial)=(( A(:,trial)-(C_.*theta_comp2_wHy) )>0 ).*( A(:,trial)-(C_.*theta_comp2_wHy));
            %                         codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
            %                        end
            %
            %                        InhAbs_CL=mean(codingLevelDummy);
            %                       depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
            %                       depsi1_dy(isnan(depsi1_dy))=0;
            %                       depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_comp2_wHy,1,odors*numtrainingSamples));
            %                       %eta=10.*mean(Y_d,2);
            %                       eta=0.5;
            %                       Grad= ((InhAbs_CL)-0.20)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
            %                       C_=C_ - (eta*Grad);
            %
            %                       if (C_<0)
            %                          error('the scale factor in the blue model wHY is -ve!!')
            %                       end
            %
            %                 % replicating the sparsity level of the KCs in real MB network
            %                 % CL=10%
            %                       eta_2=0.00000001;
            %                       for trial = 1:(odors*numtrainingSamples)
            %
            %                         Activations(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial);
            %                         Y_(:,trial)=(( Activations(:,trial)-(APLgains_wHy)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_wHy))>0 ).*( Activations(:,trial)-APLgains_wHy*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_wHy));
            %                         codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
            %
            %                       end
            %                      CL_=mean(codingLevelDummy);
            %                      dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
            %                      dsig_dy(isnan(dsig_dy))=0;
            %                      dAct_dalpha= sum(Activations,1);
            %
            %                       if ( any(isinf(dAct_dalpha)) )
            %
            %                           stopp=1;
            %                       end
            %
            %                       dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
            %                       Grad_alpha= ((CL_)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
            %                       APLgains_wHy= APLgains_wHy- eta_2*(Grad_alpha);
            %
            %
            %                 %now, third optimization step: finding the optimum weights
            %                 %to equalize the average activity levels among active KCs
            %                 %ONLY. this rule is the original one derived from the
            %                 %error function; without any heuristic for the silent KCs.
            %
            %                       for trial = 1:(odors*numtrainingSamples)
            %
            %                        Activations(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial );
            %                         Y_(:,trial)=(( Activations(:,trial)-(APLgains_wHy)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_wHy))>0 ).*( Activations(:,trial)-APLgains_wHy*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_wHy));
            %
            %                       end
            %                     Conn2=repmat(Conn,1,1,odors*numtrainingSamples);
            %                     Xik=reshape(PNtrials(:,1:odors*numtrainingSamples),m,1,odors*numtrainingSamples);
            %                     Xik_=repmat(Xik,1,n,1);
            %                     H_Y=(Y_>0);
            %                     H_Y_=repmat(H_Y,1,1,m);
            %                     H_Y_2=permute(H_Y_,[3 1 2]);
            %                     filt_=  (H_Y_2.* ((1-(APLgains_wHy)) .* (Xik_.* Conn2)) ) ;
            %                     filt_Xjm= mean(filt_,3);
            %
            %                    % change weight for already active cells,
            %                    % if a KC is silent, leave it as is. and don't change its excitability.
            %
            %                   avgAKcs=mean(Y_,2);
            %                   errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);
            %                   mask=filt_Xjm;
            %                   thisW_ActivityBasedComp_wHy= thisW_ActivityBasedComp_wHy-(0.05).*((1.*(mask.*errorInActivity)));
            %
            %
            %                   if (~isempty(find(isinf(thisW_ActivityBasedComp_wHy) )))
            %                       g=1;
            %                   end
            %                   if (any(avgAKcs==0))
            %                       silentKc=1;
            %                   end
            %
            %                   thisW_ActivityBasedComp_wHy(find(thisW_ActivityBasedComp_wHy<0))=0;
            %
            %                   % check if the constraints are satisfied:
            %
            %                   for trial = 1:(odors*numtrainingSamples)
            %
            %                     ActivationsDummy(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial);
            %                     Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_comp2_wHy))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_comp2_wHy));
            %                     codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
            %                   end
            %                  InhAbs_CL=mean(codingLevelDummy);
            %
            %                   for trial = 1:(odors*numtrainingSamples)
            %
            %                     Activations(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial);
            %                     Y_(:,trial)=(( Activations(:,trial)-(APLgains_wHy)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_wHy))>0 ).*( Activations(:,trial)-APLgains_wHy*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_wHy));
            %                     codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
            %
            %                   end
            %                  CL_=mean(codingLevelDummy);
            %                  avgAKcs=mean(Y_,2);
            %
            %                  % break off the loop only if the sparsity constraints are
            %                  % met, and all the firing KCs have the same average
            %                  % activity, i.e the KCs which won't meet target activity
            %                  % level A0 have to be the silent KCs only
            %
            %                  strayKCs= avgAKcs(find(abs(avgAKcs-A0)>epsilon));
            %                  conditions=  ( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01 ) & (all(strayKCs==0));
            %                 end
            %
            %                 theta_comp2_wHy=(C_.*theta_comp2_wHy);
            
            %% tuning theta for equalizing KCs firing probabilities: Kennedy's 2019 inspired model
            % tic
            %                 Activations =zeros(n,odors*numtrainingSamples);
            %                 ActivationsDummy= zeros(n, odors*numtrainingSamples);
            %                 APLgains(4)=0;
            %                 theta_Kenn= 1+ 1.*rand(2000,1); %% avoid negative values of theta
            %                 conditions=0;
            %                 A=zeros(n,odors*numtrainingSamples);
            %                 Y_d=zeros(n,odors*numtrainingSamples);
            %                 codingLevelDummy=[];
            %
            %                 while(~conditions)
            %
            %                 % with inhibition gain absent
            %
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     A(:,trial) = thisW_Kennedy'*PNtrials(:,trial);
            %                     Y_d(:,trial)=(( A(:,trial)-theta_Kenn)>0 ).*( A(:,trial)-theta_Kenn);
            %                 end
            %
            %                 InhAbs_CL= sum((Y_d>0),2)./(odors*numtrainingSamples);
            %                 depsi1_dy=(exp(1.*Y_d)./((1+exp(1.*Y_d)).^2));
            %                 depsi1_dy(isnan(depsi1_dy))=0;
            %                 depsi1_dtheta= -(Y_d>0).* depsi1_dy;
            %                 eta=10.*mean(Y_d,2);
            %
            %
            %                 Grad= ((InhAbs_CL)-0.20)*(1/(odors*numtrainingSamples)).*(sum(depsi1_dtheta,2));
            %                 theta_Kenn=theta_Kenn - (eta.*Grad);
            %
            %                 % stuck in a minima, low activity! then set its thetas to
            %                 % decay exponentially, by factor of 0.8.
            %                 ind_lowGrad=find((Grad)<=1e-6); %+ve, -ve and zero gradients.
            %                 pstve_zerolowGrad= ind_lowGrad( find(Grad(ind_lowGrad)>=0)); %exclude -ves only
            %                 pstve_zerolowGrad= pstve_zerolowGrad(InhAbs_CL(pstve_zerolowGrad)==0); %scale down thetas for KCs
            %                 theta_Kenn(pstve_zerolowGrad)= theta_Kenn(pstve_zerolowGrad).*0.8; % drop those thetas by 0.8
            %
            %                 % stuck in a minima, high activity! then set its thetas to
            %                 % grow exponentially by factor of 1.2.
            %
            %                 ind_lowGrad=find((Grad)<=1e-6); %+ve, -ve and zero gradients.
            %                 ngtve_zerolowGrad= ind_lowGrad( find(Grad(ind_lowGrad)<=0)); %exclude +ves
            %                 ngtve_zerolowGrad= ngtve_zerolowGrad(InhAbs_CL(ngtve_zerolowGrad)==1); %scale up thetas for these KCs, highly active ones
            %                 theta_Kenn(ngtve_zerolowGrad)= theta_Kenn(ngtve_zerolowGrad).*1.2; % scale up those thetas by 0.8
            %
            %                 % hard limit the thetas to be +ve or zero.
            %                 theta_Kenn(theta_Kenn<0)=0.01;
            %
            %
            %                 eta_2=0.0000001;
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     Activations(:,trial) = thisW_Kennedy'*PNtrials(:,trial);
            %                     Y_(:,trial)=(( Activations(:,trial)-(APLgains(4))*repmat(sum(Activations(:,trial),1),n,1)-theta_Kenn)>0 ).*( Activations(:,trial)-APLgains(4)*repmat(sum(Activations(:,trial),1),n,1)-theta_Kenn);
            %                     codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
            %
            %                 end
            %                 CL_=mean(codingLevelDummy);
            %                 dsig_dy=(exp(1.*Y_)./((1+exp(1.*Y_)).^2));
            %                 dsig_dy(isnan(dsig_dy))=0;
            %                 dAct_dalpha= sum(Activations,1);
            %
            %                 if ( any(isinf(dAct_dalpha)) )
            %                   stopp=1;
            %                 end
            %
            %                 dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
            %                 Grad_alpha= ((CL_)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
            %                 APLgains(4)= APLgains(4)- eta_2*(Grad_alpha);
            %
            %
            %                 %% check the constraints
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     ActivationsDummy(:,trial) = thisW_Kennedy'*PNtrials(:,trial);
            %                     Y_d(:,trial)=(( ActivationsDummy(:,trial)-theta_Kenn)>0 ).*( ActivationsDummy(:,trial)-theta_Kenn);
            %                 end
            %                 InhAbs_CL= sum((Y_d>0),2)./(odors*numtrainingSamples);
            %
            %                 for trial = 1:(odors*numtrainingSamples)
            %
            %                     Activations(:,trial) = thisW_Kennedy'*PNtrials(:,trial);
            %                     Y_(:,trial)=(( Activations(:,trial)-(APLgains(4))*repmat(sum(Activations(:,trial),1),n,1)-theta_Kenn)>0 ).*( Activations(:,trial)-APLgains(4)*repmat(sum(Activations(:,trial),1),n,1)-theta_Kenn);
            %                     codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
            %
            %                 end
            %                 CL_=mean(codingLevelDummy);
            %                 conditions= all(  round(abs( round((InhAbs_CL./CL_),1) - (2.0.*ones(2000,1)) ),1) <= 0.2.*ones(2000,1) ) &( (abs(CL_-0.10)) <=0.01 );
            %
            %                 end
            toc
            
            
            %% theta tuned for activity equalization, our main model for tuning theta in Fig5
            tic
            Activations =zeros(n,odors*numtrainingSamples);
            ActivationsDummy= zeros(n, odors*numtrainingSamples);
            
            % scale the weights matrix by 2?
            thisW_Kennedy= 5.*thisW_Kennedy;
            
            APLgains(5)=0;
            theta_Activity_homeo=  1+ 1.*rand(2000,1); %% avoid negative values of theta
            conditions=0; %%
            
            A=zeros(n,odors*numtrainingSamples);
            Y_d=zeros(n,odors*numtrainingSamples);
            codingLevelDummy=[];
            A0=(0.51).*ones(n,1);
            epsilon= A0(1)*0.06;
            C_=1;
            
            t=1;
            while(~conditions)
                
                % with inhibition gain absent, make sure that
                % CL=20% with the random values for theta.
                % adjust C_theta
                A = multiprod((thisW_Kennedy)',PNtrials(:,:,1:numtrainingSamples));
                A=reshape(A,n,odors*numtrainingSamples);
                Y_d=(( A- (C_.*theta_Activity_homeo))>0 ).*( A-(C_.*theta_Activity_homeo));
                codingLevelDummy=  (sum(Y_d>0,1)/n);
                InhAbs_CL=mean(codingLevelDummy);
                depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_Activity_homeo,1,odors*numtrainingSamples));
                eta=1;
                Grad= ((InhAbs_CL)-0.20)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
                C_=C_- (eta*mean(Y_d(:))*Grad); % where is the Ybar-A0 factor?
                if (C_<0)
                    error('the scale factor in our magenta model is -ve')
                end
                
                % CL=10% constraint, adjust APL gain
                eta_2=0.000000005;
                Activations= multiprod(thisW_Kennedy',PNtrials(:,:,1:numtrainingSamples));
                Activations=reshape(Activations,n,odors*numtrainingSamples);
                Y_=(( Activations- ((repmat(APLgains(5),1,odors*numtrainingSamples)).*repmat(sum(Activations,1),n,1)) -(C_.*theta_Activity_homeo))>0 ).*( Activations-((repmat(APLgains(5),1,odors*numtrainingSamples)).*repmat(sum(Activations,1),n,1)) -(C_.*theta_Activity_homeo));
                codingLevelDummy=  (sum(Y_>0,1)/n);
                CL_=mean(codingLevelDummy);
                dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(Activations,1);
                if ( any(isinf(dAct_dalpha)) )
                    stopp=1;
                end
                dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CL_)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgains(5)= APLgains(5)- eta_2*(Grad_alpha);
                
                % now tune theta for equal average activity level.
                for trial = 1:(odors*numtrainingSamples)
                    Activations(:,trial) = (thisW_Kennedy)'*PNtrials(:,trial );
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(5))*repmat(sum(Activations(:,trial),1),n,1)- (C_.*theta_Activity_homeo))>0 ) ...
                        .* (       Activations(:,trial)- APLgains(5) *repmat(sum(Activations(:,trial),1),n,1)- (C_.*theta_Activity_homeo));
                end
                avgAKcs=mean(Y_,2);
                errorInActivity=(avgAKcs-A0);
                theta_Activity_homeo= theta_Activity_homeo-(0.01).*((-1.*C_.*(errorInActivity)));
                theta_Activity_homeo(theta_Activity_homeo<0)=0;
                
                % check the constraints
                for trial = 1:(odors*numtrainingSamples)
                    ActivationsDummy(:,trial) = (thisW_Kennedy)'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)- (C_.*theta_Activity_homeo) )>0 ).*( ActivationsDummy(:,trial)- (C_.*theta_Activity_homeo));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    Activations(:,trial) = (thisW_Kennedy)'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains(5))*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_Activity_homeo))>0 ).*( Activations(:,trial)-APLgains(5)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_Activity_homeo));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                end
                CL_=mean(codingLevelDummy);
                avgAKcs=mean(Y_,2);
                
                conditions= all(abs(avgAKcs-A0)<epsilon) &( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01);
                
                InhAbs_CL_vec(t)=InhAbs_CL;
                avgAct(t,:)= (avgAKcs);
                CL_vec(t)= CL_;
                if mod(t,20)==0
                    disp(length(find(abs(avgAKcs-A0)>epsilon)))
                    disp(CL_)
                end
                t=t+1;
                
            end
            disp('time to run magenta model:');
            toc
            theta_Activity_homeo=C_.*theta_Activity_homeo;
            
            %% inhibitory plasticity compensation
            tic
            APLtrajectory=zeros(2000,1);
            CLtrajectory=[];
            AvgAKCstrajectory=zeros(2000,1);
            thisW_ActivityBasedComp_inhibitionPlast= 5.*thisW;
            C_=1;
            cw=1;
            Activations =zeros(n,odors*numtrainingSamples);
            ActivationsDummy= zeros(n, odors*numtrainingSamples);
            APLgains_model6=zeros(2000,1);
            T=5;
            A0=(0.51).*ones(n,1);
            epsilon= A0(1)*0.06;
            
            theta_inhibitionPlast= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            theta_inhibitionPlast(theta_inhibitionPlast>70)=70;
            theta_inhibitionPlast(theta_inhibitionPlast<0.01)=0.01;
            
            t=1;
            conditions=0;
            
            A=zeros(n,odors*numtrainingSamples);
            Y_d=zeros(n,odors*numtrainingSamples);
            codingLevelDummy=[];
            
            while(~conditions)
                
                
                A = multiprod( ( thisW_ActivityBasedComp_inhibitionPlast)',PNtrials(:,:,1:numtrainingSamples));
                A=reshape(A,n,odors*numtrainingSamples);
                
                Y_d=(( A-(C_.*theta_inhibitionPlast))>0 ).*( A-(C_.*theta_inhibitionPlast));
                codingLevelDummy=  (sum(Y_d>0,1)/n);
                
                InhAbs_CL=mean(codingLevelDummy);
                
                depsi1_dy=(exp(1.*Y_d)./((1+exp(1.*Y_d)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                
                depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_inhibitionPlast,1,odors*numtrainingSamples));
                
                
                Grad= (((InhAbs_CL)-0.19)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) )));
                
                C_=C_ - (1*mean(Y_d(:)).*(Grad));
                
                
                if (C_<0 )
                    error('the scale factor in green model is -ve')
                end
                
                Activations = multiprod( ( thisW_ActivityBasedComp_inhibitionPlast)',PNtrials(:,:,1:numtrainingSamples));
                Activations= reshape(Activations,n,odors*numtrainingSamples);
                
                Y_= (Activations-(repmat(APLgains_model6,1,odors*numtrainingSamples).*(repmat(sum(Activations),n,1) ))-(C_.*theta_inhibitionPlast)>0 ).*( Activations-(repmat(APLgains_model6,1,odors*numtrainingSamples).*(repmat(sum(Activations),n,1)))-(C_.*theta_inhibitionPlast));
                
                codingLevelDummy=  (sum(Y_>0,1)/n);
                CL_=mean(codingLevelDummy);
                avgAKcs=mean(Y_,2);
                errorInActivity=avgAKcs-A0;
                dsig_dy=(exp(1.*Y_)./((1+exp(1.*Y_)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                xk_=repmat(sum(Activations),n,1);
                D_yj_alphaj= (-1.*mean(xk_,2));
                dAct_dalpha= sum(Activations,1);
                dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                dYik_dalphai= -(repmat(dAct_dalpha,n,1));
                
                
                eta_o1= 0.02;
                eta_o2= 0.0000001;
                
                % for some KCs change in alpha for activity equalization can cancel out
                % this happen in few KCs towards approaching the global minima for
                % all KCs, the to avoid stucking in a local minima, we add a small noise
                % on top of the calculated gradient.
                
                Grad_alpha1= ( ((eta_o1)) .*((CL_)-0.10)*(1/(n*odors*numtrainingSamples)).*(sum(dsig_dalpha,2)) ) + (0.0001*randn(n,1));
                
                Grad_alpha2=( eta_o2.* (errorInActivity).*(D_yj_alphaj) ) + (0.0001*randn(n,1));
                
                Grad_alpha= Grad_alpha1+Grad_alpha2 ;
                
                APLgains_model6= APLgains_model6- 0.001.*((Grad_alpha));
                
                
                %% check constraints
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    ActivationsDummy(:,trial) = (thisW_ActivityBasedComp_inhibitionPlast)'*PNtrials(:,trial);
                    Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_inhibitionPlast))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_inhibitionPlast));
                    codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
                end
                
                InhAbs_CL=mean(codingLevelDummy);
                
                for trial = 1:(odors*numtrainingSamples)
                    
                    Activations(:,trial) = (thisW_ActivityBasedComp_inhibitionPlast)'*PNtrials(:,trial);
                    Y_(:,trial)=(( Activations(:,trial)-(APLgains_model6.*sum(Activations(:,trial),1))-(C_.*theta_inhibitionPlast))>0 ).*( Activations(:,trial)-(APLgains_model6.*sum(Activations(:,trial),1))-(C_.*theta_inhibitionPlast));
                    codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
                    
                end
                CL_=mean(codingLevelDummy);
                
                avgAKcs=mean(Y_,2);
                
                conditions= all(abs(avgAKcs-A0)<epsilon)  &( (abs(round(CL_,3)-0.10)) <=0.015) & ( abs( (InhAbs_CL/CL_) - 2.0) <0.3 );
                
                InhAbs_CL_vec(t)=InhAbs_CL;
                avgAct(t,:)= (avgAKcs);
                CL_vec(t)= CL_;
                if mod(t,20)==0
                    disp(length(find(abs(avgAKcs-A0)>epsilon)))
                    disp(CL_)
                end
                t=t+1;
                
            end
            theta_inhibitionPlast=(C_.*theta_inhibitionPlast);
            thisW_ActivityBasedComp_inhibitionPlast= (cw.*thisW_ActivityBasedComp_inhibitionPlast);
            
            disp('time to run green model:');
            toc
            
            %                save the parameters for each fly, connectivity weights,
            %                spiking thresholds and inhibition Gains in all the models
            
            save( strcat(CurrDir,'/TunedFlies_allModels_HOInp_wRound',['_fly_wNoise',num2str(randomTrials),num2str(noiseScale)]) ,'APLgains',...
                'PNtrials','thisW_HomogModel','thetaH_Ftheta', ...
                'thisW','thetaS', 'thisW_equalizedModel', 'InhibitionGain','theta',...
                'APLgains_model6','thisW_ActivityBasedComp_inhibitionPlast','theta_inhibitionPlast',...
                'thisW_Kennedy','theta_Activity_homeo','classAction1','APLgains_noxjk','thisW_ActivityBasedComp_noxjk' ,'theta_comp2_noxjk',...
                'PnToKc');
            
            
            Activations=zeros(n,odors*numTrials);
            ActivationsEqualized=zeros(n,odors*numTrials);
            ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numTrials);
            Activations_comp2= zeros(n,odors*numTrials);
            Activations_comp2_wHy= zeros(n,odors*numTrials);
            Activations_comp2_noxjk=zeros(n,odors*numTrials);
            Activations_Kenn= zeros(n,odors*numTrials);
            Activations_theta_activity_hoemeo= zeros(n,odors*numTrials);
            Activations_inhibPlast=zeros(n,odors*numTrials);
            
            Y=zeros(n,odors*numTrials);
            YEqualized=zeros(n,odors*numTrials);
            YHomogdummy_Ftheta=zeros(n,odors*numTrials);
            Y_comp2= zeros(n,odors*numTrials);
            Y_Kenn=zeros(n,odors*numTrials);
            Y_inhibPlast= zeros(n,odors*numTrials);
            Y_theta_activity_homeo=zeros(n,odors*numTrials);
            Y_comp2_noxjk=zeros(n,odors*numTrials);
            Y_comp2_wHy= zeros(n,odors*numTrials);
            
            
            for trial = 1:(odors*numTrials)
                
                
                ActivationsHomogenousdummy_Ftheta(:,trial) = thisW_HomogModel'*PNtrials(:,trial  );
                YHomogdummy_Ftheta(:,trial)=(( ActivationsHomogenousdummy_Ftheta(:,trial)-(APLgains(2))*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomogenousdummy_Ftheta(:,trial)-APLgains(2)*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);
                
                Activations_theta_activity_hoemeo(:,trial) = thisW_Kennedy'*PNtrials(:,trial );
                Y_theta_activity_homeo(:,trial)=(( Activations_theta_activity_hoemeo(:,trial)-(APLgains(5) )*repmat(sum(Activations_theta_activity_hoemeo(:,trial),1),n,1)-theta_Activity_homeo)>0 ).*( Activations_theta_activity_hoemeo(:,trial)-APLgains(5)*repmat(sum(Activations_theta_activity_hoemeo(:,trial),1),n,1)-theta_Activity_homeo);
                
                Activations(:,trial) = thisW'*PNtrials(:,trial );
                Y(:,trial)=(( Activations(:,trial)-(APLgains(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgains(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
                
                ActivationsEqualized(:,trial) = thisW_equalizedModel'*PNtrials(:,trial );
                YEqualized(:,trial)=(( ActivationsEqualized(:,trial)-(InhibitionGain)*repmat(sum(ActivationsEqualized(:,trial),1),n,1)-theta)>0 ).*( ActivationsEqualized(:,trial)-InhibitionGain*repmat(sum(ActivationsEqualized(:,trial),1),n,1)-theta);
                
                %                     Activations_comp2(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial );
                %                     Y_comp2(:,trial)=(( Activations_comp2(:,trial)-(APLgains(3) )*repmat(sum(Activations_comp2(:,trial),1),n,1)-theta_comp2)>0 ).*( Activations_comp2(:,trial)-APLgains(3)*repmat(sum(Activations_comp2(:,trial),1),n,1)-theta_comp2);
                
                Activations_comp2_noxjk(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial );
                Y_comp2_noxjk(:,trial)=(( Activations_comp2_noxjk(:,trial)-(APLgains_noxjk )*repmat(sum(Activations_comp2_noxjk(:,trial),1),n,1)-theta_comp2_noxjk)>0 ).*( Activations_comp2_noxjk(:,trial)-APLgains_noxjk*repmat(sum(Activations_comp2_noxjk(:,trial),1),n,1)-theta_comp2_noxjk);
                
                %                     Activations_comp2_wHy(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial );
                %                     Y_comp2_wHy(:,trial)=(( Activations_comp2_wHy(:,trial)-(APLgains_wHy )*repmat(sum(Activations_comp2_wHy(:,trial),1),n,1)-theta_comp2_wHy)>0 ).*( Activations_comp2_wHy(:,trial)-APLgains_wHy*repmat(sum(Activations_comp2_wHy(:,trial),1),n,1)-theta_comp2_wHy);
                %
                %                     Activations_Kenn(:,trial) = thisW'*PNtrials(:,trial );
                %                     Y_Kenn(:,trial)=(( Activations_Kenn(:,trial)-(APLgains(4) )*repmat(sum(Activations_Kenn(:,trial),1),n,1)-theta_Kenn)>0 ).*( Activations_Kenn(:,trial)-APLgains(4)*repmat(sum(Activations_Kenn(:,trial),1),n,1)-theta_Kenn);
                
                Activations_inhibPlast(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials(:,trial );
                Y_inhibPlast(:,trial)=(( Activations_inhibPlast(:,trial)-(APLgains_model6 ).*repmat(sum(Activations_inhibPlast(:,trial),1),n,1)-theta_inhibitionPlast)>0 ).*( Activations_inhibPlast(:,trial)-APLgains_model6.*repmat(sum(Activations_inhibPlast(:,trial),1),n,1)-theta_inhibitionPlast);
                
            end
            
            
            for l_r=1:lrs
                
                WopAllOdours=1*rand(n,2);
                WopAllOdoursEqualized= WopAllOdours;
                
                WopAllOdoursHomog_Ftheta=WopAllOdours;
                WopAllOdoursInhPlast=WopAllOdours;
                
                WopAllOdoursEqualizedComp2=WopAllOdours;
                
                WopFromPNs= 1*rand(m,2);
                
                WopAllOdoursKenn=WopAllOdours;
                
                WopAllOdoursThetaHomeo= WopAllOdours;
                
                WopAllOdoursEqualizedComp2_noxjk=WopAllOdours;
                WopAllOdoursEqualizedComp2_wHy= WopAllOdours;
                
                
                alpha=0.000001* (10^((l_r)));
                
                c=1;
                ceq=1;
                ch=1;
                
                
                YHomog_Fthetatemp=reshape(YHomogdummy_Ftheta,n,odors,numTrials);
                
                Ytemp= reshape(Y,n,odors,numTrials);
                
                YEqualizedtemp=reshape(YEqualized,n,odors,numTrials);
                
                Y_comp2temp= reshape(Y_comp2,n,odors,numTrials);
                
                Y_Kenntemp= reshape(Y_Kenn,n,odors,numTrials);
                
                Y_theta_activity_homeotemp= reshape(Y_theta_activity_homeo,n,odors,numTrials);
                
                Y_inhibPlasttemp= reshape(Y_inhibPlast,n,odors,numTrials);
                
                Y_comp2_noxjktemp= reshape(Y_comp2_noxjk,n,odors,numTrials);
                
                Y_comp2_wHytemp= reshape(Y_comp2_wHy,n,odors,numTrials);
                
                
                
                % rescale the KC responses to be from 0-1
                Ytemp=rescale(Ytemp);
                YEqualizedtemp=rescale(YEqualizedtemp);
                YHomog_Fthetatemp=rescale(YHomog_Fthetatemp);
                Y_comp2temp=rescale(Y_comp2temp);
                Y_Kenntemp= rescale(Y_Kenntemp);
                Y_inhibPlasttemp=rescale(Y_inhibPlasttemp);
                Y_theta_activity_homeotemp=rescale(Y_theta_activity_homeotemp);
                Y_comp2_noxjktemp= rescale(Y_comp2_noxjktemp);
                Y_comp2_wHytemp=rescale(Y_comp2_wHytemp);
                
                YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
                
                Ytr=Ytemp(:,:,1:numtrainingSamples);
                
                YEqualizedtr=YEqualizedtemp(:,:,1:numtrainingSamples);
                
                Y_comp2tr= Y_comp2temp(:,:,1:numtrainingSamples);
                
                Y_Kenntr= Y_Kenntemp(:,:,1:numtrainingSamples);
                
                Y_inhibPlasttr= Y_inhibPlasttemp(:,:,1:numtrainingSamples);
                
                Y_theta_activity_homeotr= Y_theta_activity_homeotemp(:,:,1:numtrainingSamples);
                
                Y_comp2_noxjktr= Y_comp2_noxjktemp(:,:,1:numtrainingSamples);
                
                Y_comp2_wHytr= Y_comp2_wHytemp(:,:,1:numtrainingSamples);
                
                
                %% learning from a preceptron in the output layer
                
                
                
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
                        
                        delta =  exp( -(alpha/mean(Y_comp2_noxjktr(:)))* sum(Y_comp2_noxjktr(:,odour,:),3) );
                        
                        WopAllOdoursEqualizedComp2_noxjk(:,2)= WopAllOdoursEqualizedComp2_noxjk(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_comp2_noxjktr(:)))* sum(Y_comp2_noxjktr(:,odour,:),3) );
                        WopAllOdoursEqualizedComp2_noxjk(:,1)= WopAllOdoursEqualizedComp2_noxjk(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_comp2_wHytr(:)))* sum(Y_comp2_wHytr(:,odour,:),3) );
                        
                        WopAllOdoursEqualizedComp2_wHy(:,2)= WopAllOdoursEqualizedComp2_wHy(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_comp2_wHytr(:)))* sum(Y_comp2_wHytr(:,odour,:),3) );
                        WopAllOdoursEqualizedComp2_wHy(:,1)= WopAllOdoursEqualizedComp2_wHy(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_Kenntr(:)))* sum(Y_Kenntr(:,odour,:),3) );
                        
                        WopAllOdoursKenn(:,2)= WopAllOdoursKenn(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_Kenntr(:)))* sum(Y_Kenntr(:,odour,:),3) );
                        WopAllOdoursKenn(:,1)= WopAllOdoursKenn(:,1) .*delta;
                        
                    end
                    
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_inhibPlasttr(:)))* sum(Y_inhibPlasttr(:,odour,:),3) );
                        
                        WopAllOdoursInhPlast(:,2)= WopAllOdoursInhPlast(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_inhibPlasttr(:)))* sum(Y_inhibPlasttr(:,odour,:),3) );
                        WopAllOdoursInhPlast(:,1)= WopAllOdoursInhPlast(:,1) .*delta;
                        
                    end
                    
                    if( ~ isempty(find(classAction1==odour)) )
                        
                        delta =  exp( -(alpha/mean(Y_theta_activity_homeotr(:)))* sum(Y_theta_activity_homeotr(:,odour,:),3) );
                        
                        WopAllOdoursThetaHomeo(:,2)= WopAllOdoursThetaHomeo(:,2) .*delta;
                        
                    else
                        
                        delta = exp(- (alpha/mean(Y_theta_activity_homeotr(:)))* sum(Y_theta_activity_homeotr(:,odour,:),3) );
                        WopAllOdoursThetaHomeo(:,1)= WopAllOdoursThetaHomeo(:,1) .*delta;
                        
                    end
                    
                    
                end
                
                %% performance as a function of the strictness of the decision making
                % this strictness is dictated by C in the soft-max function.
                % so given the same fly, same task, and after learning measure the performance as f(c)
                
                for c=1:Crange
                    
                    
                    C=C_SoftMax*(10^c);
                    
                    [acc,accEq,accEq2_noxjk,accKenn,accInhPlast,accThetaHomeo]=testingModels_accuracies(C,WopAllOdours,WopAllOdoursEqualized, WopAllOdoursEqualizedComp2_noxjk, WopAllOdoursKenn, WopAllOdoursInhPlast,...
                        WopAllOdoursThetaHomeo, classAction1,numTrials,numtrainingSamples,Ytemp,YEqualizedtemp,Y_comp2_noxjktemp,Y_Kenntemp,Y_inhibPlasttemp,Y_theta_activity_homeotemp);
                    
                    test_p_raEq(randomTrials, l_r,c)=accEq;
                    test_p_ra(randomTrials, l_r,c)=acc;
                    %test_p_raEq2(randomTrials, l_r,c)=accEq2;
                    test_p_raKenn(randomTrials, l_r,c)= accKenn;
                    test_p_raInhPlast(randomTrials, l_r,c)= accInhPlast;
                    test_p_ra_thetaHomeo(randomTrials, l_r,c)= accThetaHomeo;
                    test_p_raEq2_noxjk(randomTrials, l_r,c)= accEq2_noxjk;
                    %test_p_raEq2_wHy(randomTrials, l_r,c)= accEq2_wHy;
                    
                    [accH1]=HomogenousModel_KernelTesting(C,WopAllOdoursHomog_Ftheta,PNtrials, thetaH_Ftheta,InhibitionGain, APLgains,classAction1,numTrials,numtrainingSamples,YHomog_Fthetatemp);
                    
                    test_p_raH_FixedTheta(randomTrials, l_r,c)=accH1;
                    
                    
                    
                    
                end
                
            end
            
        end
        
        
    end
% end

%% save the following in the MATLAB current directory.

%% save the models' accuracy scores:
% save('test_p_ra_withround_HOinp.mat','test_p_ra');
% save('test_p_raEq_withround_HOinp.mat','test_p_raEq');
% save('test_p_raH_FixedTheta_withround_HOinp.mat','test_p_raH_FixedTheta');
% save('test_p_raInhPlast_withround_HOinp.mat','test_p_raInhPlast');
% save('test_p_ra_thetaHomeo_withround_HOinp.mat','test_p_ra_thetaHomeo');
% save('test_p_raEq2_noxjk_withround_HOinp.mat','test_p_raEq2_noxjk');


%% End of the code---------------------------------------------------------

delete(gcp('nocreate'));



