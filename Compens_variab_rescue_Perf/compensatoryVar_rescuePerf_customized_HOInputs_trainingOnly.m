% allows running the training only of the model flies saved from
% compensatoryVar_rescuePerf_customized_HOInputs

%% Start of the code

clc
clear all;

CurrDir= pwd;

load('hallem_olsen.mat');
load('PW_given_N.mat');
load('PW_given_theta_and_n.mat');
load('W_PN_KC.mat');
load ('P_n.mat');
ll=20;
numtrainingSamples=15; %% artificial plus linearly dependent responses
modss=1;
mulOd=110;
lrs=7;
logLRs = [-5 -4 -3 -2.75 -2.5 -2.25 -2 -1 0 1];
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

% % start a local parallel pool of 3 threads to speed up simulation time
% p=parpool(3);
% p.IdleTimeout = 1000;
% parfevalOnAll(@maxNumCompThreads,0,3)


for mods=1:modss
    
    odors= mulOd*(mods/modss);
    x=PN;
    
    for randomTrials=1:ll
        
        n =2000; % number of neurons in the hidden layer
        m=24;  %number of dimensions in the input data
        
        randomTrials
        
        for noiseScale=1:NScales
            
            
            
            %                save the parameters for each fly, connectivity weights,
            %                spiking thresholds and inhibition Gains in all the models
            
            load( strcat(CurrDir,'/TunedFlies_allModels_HOInp_wRound',['_fly_wNoise',num2str(randomTrials),num2str(noiseScale)]));
            
            
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
            
            
            for l_r=1:length(logLRs)
                
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
                
                
                alpha=10^logLRs(l_r);%0.000001* (10^((l_r)));
                
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
                
                %% perfromance as a function of the strictness of the decision making
                %% this strictness is dictated by C in the soft-max function.
                %% so given the same fly, same task, and after learning measure the performance as f(c)
                
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
end

%% save the following in the MATLAB current directory.
%cd pwd
%% save the models' accuracy scores:
save('test_p_ra_withround_HOinp.mat','test_p_ra');
save('test_p_raEq_withround_HOinp.mat','test_p_raEq');
save('test_p_raH_FixedTheta_withround_HOinp.mat','test_p_raH_FixedTheta');
save('test_p_raInhPlast_withround_HOinp.mat','test_p_raInhPlast');
save('test_p_ra_thetaHomeo_withround_HOinp.mat','test_p_ra_thetaHomeo');
save('test_p_raEq2_noxjk_withround_HOinp.mat','test_p_raEq2_noxjk');


%% End of the code---------------------------------------------------------

delete(gcp('nocreate'));




