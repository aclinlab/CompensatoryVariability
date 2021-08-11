%% Date: 11 August 2021
% this code reproduces the results in 3F
% S2(H,I)
% this code highlights the effect of the sparsest
% KCs in rescuing the memory performance at dense coding levels.
%%


clear all
clc
load('hallem_olsen.mat');
numtrainingSamples=15; %% artificial plus linearly dependent responses 
mulOd=20;
% mulOd=100;
lrs=10;
Crange=1;           
odors=mulOd;
numTrials = 30;
n =2000; % number of neurons in the hidden layer
m=24; 
C_SoftMax=1;
NScales=1;
LRs=10.^[-5 -4 -3 -2.75 -2.5 -2.25 -2 -1 0 1]; 

Spar_rand_ind=[];
Spar_homog_ind=[];
Spec_rand_ind=[];
Spec_homog_ind=[];

CL=0.9;

for randomTrials=1:50
    
    % load the tuned fly models, tuned on the subset of 20 odors task at
    % CL=0.9
    load(strcat(['CL_performance_random_homog20_',num2str(CL),'_',num2str(randomTrials),'.mat']));

    % load the tuned fly models, tuned on the 100 odors task at CL=0.9
    % load(strcat(['CL_performance_random_homog100_',num2str(CL),'_',num2str(randomTrials),'.mat']));
    
    
    Activations=zeros(n,odors*numTrials);
    ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numTrials);
    Y=zeros(n,odors*numTrials);
    YHomog_Ftheta=zeros(n,odors*numTrials);
    
    for trial = 1:(odors*numTrials)
       ActivationsHomogenousdummy_Ftheta(:,trial) = thisW_HomogModel'*PNtrials(:,trial  );
       YHomog_Ftheta(:,trial)=(( ActivationsHomogenousdummy_Ftheta(:,trial)-(APLgainP(2))*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomogenousdummy_Ftheta(:,trial)-APLgainP(2)*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);                   

       Activations(:,trial) = thisW'*PNtrials(:,trial );
       Y(:,trial)=(( Activations(:,trial)-(APLgainP(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgainP(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
    end

        %% calculate LT sparsity:
         
        Kpats=odors*numTrials;

        Spar=(1/(1-(1/Kpats)))*(1-(( sum( (Y./Kpats),2 ) ).^2./(sum( ( (Y.^2)./Kpats),2 ) )));

        HomogenousSpar= (1/(1-(1/Kpats)))*(1-((sum( (YHomog_Ftheta./Kpats),2)).^2./(sum( ( (YHomog_Ftheta.^2)./Kpats),2 ))));

        lifeTime_Spar(:,randomTrials)= Spar;
        
        H_lifeTime_Spar(:,randomTrials)= HomogenousSpar;

        % sort the sparsity from smallest to largest values.
        [v1,i1_sparsest]=sort(Spar);
        [v2,i2_sparsest]=sort(HomogenousSpar);
        Y_rem_sparseKCs=Y;
        YHomog_rem_sparseKCs=YHomog_Ftheta;
        
        % filter the sorted indcies out of the nan KCs, silent KCs:
        i1_sparsest(isnan(Spar(i1_sparsest)))=[];
        i2_sparsest(isnan(HomogenousSpar(i2_sparsest)))=[];

        %% save the IDs of the sparsest Kcs indcies per fly
           Spar_rand_ind(:,randomTrials)= i1_sparsest(end-199:end);
           Spar_homog_ind(:,randomTrials)= i2_sparsest(end-199:end);
        %%
        
        % replace the top 10% sparse KCs with useless KCs;'silent' and
        % 'always active' neurons.
        % while doing that, we make sure that the CL remains 0.9 after 
        % the KCs responses matrix is amputated with useless responses:
        
        Y_rem_sparseKCs(i1_sparsest(end-199:end),:)=[];
        YHomog_rem_sparseKCs(i2_sparsest(end-199:end),:)=[];
        
        % CL of the networks after removing the top 10% sparsest KCs.
        CL_rem_sparseKCs= mean(sum(Y_rem_sparseKCs>0)./1800);
        CL_Homog_rem_sparseKCs=mean(sum(YHomog_rem_sparseKCs>0)./1800);
        
        % amputation with 'always active' KCs:
        % number of useless always active KCs= 9(1-CL)*200
        
        KCs_alwaysActiv=round((9*(1-CL_rem_sparseKCs))*200);
        Homog_KCs_alwaysActiv=round((9*(1-CL_Homog_rem_sparseKCs))*200);
        ValueMax=max(Y(:));
        ValueMax_H=max(YHomog_Ftheta(:));
        Y(i1_sparsest(end-199:end),:)=[];
        YHomog_Ftheta(i2_sparsest(end-199:end),:)=[];
        Y(end+1:end+KCs_alwaysActiv,:)=repmat(ValueMax.*ones(1,odors*numTrials),KCs_alwaysActiv,1);
        YHomog_Ftheta(end+1:end+Homog_KCs_alwaysActiv,:)=repmat(ValueMax_H.*ones(1,odors*numTrials),Homog_KCs_alwaysActiv,1);
       
        % amputation with 'silent' KCs: 
        KCs_silent=200-(KCs_alwaysActiv);
        Homog_KCs_silent=200-Homog_KCs_alwaysActiv;
        Y(end+1:end+KCs_silent,:)=repmat(zeros(1,odors*numTrials),KCs_silent,1);
        YHomog_Ftheta(end+1:end+Homog_KCs_silent,:)=repmat(zeros(1,odors*numTrials),Homog_KCs_silent,1);
        
        % recalculate LT sparsity of the amputated networks
        Kpats=odors*numTrials;
        Spar=(1/(1-(1/Kpats)))*(1-(( sum( (Y./Kpats),2 ) ).^2./(sum( ( (Y.^2)./Kpats),2 ) )));
        HomogenousSpar= (1/(1-(1/Kpats)))*(1-((sum( (YHomog_Ftheta./Kpats),2)).^2./(sum( ( (YHomog_Ftheta.^2)./Kpats),2 ))));
        lifeTime_Spar(:,randomTrials)= Spar;        
        H_lifeTime_Spar(:,randomTrials)= HomogenousSpar;
        H_LTSpar(randomTrials)=std(HomogenousSpar(~isnan(HomogenousSpar)));
        S_LTSpar(randomTrials)=std(Spar(~isnan(Spar)));
       
         
%%          training and testing the amputated networks:
        
        for l_r=1:lrs  
               
            WopAllOdours=1*rand(n,2);              
            WopAllOdoursHomog_Ftheta=WopAllOdours;
            alpha=LRs(l_r);
            c=1;
            ceq=1;
            ch=1;
            YHomog_Fthetatemp=reshape(YHomog_Ftheta,n,odors,numTrials);
            Ytemp= reshape(Y,n,odors,numTrials);

            % rescale the KC responses to be from 0-1
            Ytemp=rescale(Ytemp);
            YHomog_Fthetatemp=rescale(YHomog_Fthetatemp);
            YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
            Ytr=Ytemp(:,:,1:numtrainingSamples);

                
            % learning from a preceptron in the output layer

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

             end
                
                
             for c=1:Crange
                 
                
                 C=C_SoftMax*(10^c);
               
                 [acc]=testingRandomModel_accuracy(C,WopAllOdours,classAction1,numTrials,numtrainingSamples,Ytemp);
                
                test_p_ra(randomTrials, l_r)=acc;
                
                
              [accH1]=HomogenousModel_KernelTesting(C,WopAllOdoursHomog_Ftheta,PNtrials, thetaH_Ftheta,classAction1,numTrials,numtrainingSamples,YHomog_Fthetatemp);

                test_p_raH_FixedTheta(randomTrials, l_r)=accH1;  
                
             end
      
        end
       
             % get the ang. distance after discarding the most sparse KCs            
             S_PW_AngDist(:,:,randomTrials)=AngDist_Pw_Calc(Y);
             H_PW_AngDist(:,:,randomTrials)=AngDist_Pw_Calc(YHomog_Ftheta);
       
end

% save the models performances & angular distances
save( strcat(['Excl_mostSparseKCs_rand_homogModels_',num2str(odors),'odors.mat']),'test_p_ra','test_p_raH_FixedTheta','S_PW_AngDist','H_PW_AngDist');

