%% Date: 11 August 2021
% this code reproduces the results in 3F
% S2(H,I)
% this code highlights the effect of the top specialised
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
    
       %% remove the top 5% specialised KCs
       
       odo=ones(1,mulOd);
       odo(classAction1)=0;
       class2=find(odo);
       
       Ytemp1=reshape(Y,2000,odors,numTrials);
       YHomog_Fthetatemp1=reshape(YHomog_Ftheta,2000,odors,numTrials);
       
       % calculate KCs valence specificity:
       V_spec_R= abs( sum(Ytemp1(:,classAction1,:),[2 3])- sum(Ytemp1(:,class2,:), [2 3]) )./ sum(Ytemp1,[2 3]); 
       V_spec_Homog= abs( sum(YHomog_Fthetatemp1(:,classAction1,:),[2 3])- sum(YHomog_Fthetatemp1(:,class2,:), [2 3]) )./ sum(YHomog_Fthetatemp1,[2 3]); 

       % sort the specificity values from smallest to largest values.        
        [v1,i1]=sort(V_spec_R);
        [v2,i2]=sort(V_spec_Homog);
       
        Y_rem_sparseKCs=Y;
        YHomog_rem_sparseKCs=YHomog_Ftheta;
        
        % filter the sorted indcies out of the nan KCs, silent KCs:
        i1(isnan(V_spec_R(i1)))=[];
        i2(isnan(V_spec_Homog(i2)))=[];
        
        %% save the IDs of the top 5% specialised Kcs indcies per fly:
        Spec_rand_ind(:,randomTrials)= i1(end-99:end);
        Spec_homog_ind(:,randomTrials)= i2(end-99:end);  
        %%
        
        % replace the top 5% specific KCs with useless KCs;'silent' and
        % 'always active' neurons.
        % while doing that, we make sure that the CL remains 0.9 after 
        % the KCs responses matrix is amputated with useless responses:
        
        Y_rem_sparseKCs(i1(end-99:end),:)=[];
        YHomog_rem_sparseKCs(i2(end-99:end),:)=[];
        
        % CL of the networks after removing the most specific KCs.
        CL_rem_sparseKCs= mean(sum(Y_rem_sparseKCs>0)./1900);
        CL_Homog_rem_sparseKCs=mean(sum(YHomog_rem_sparseKCs>0)./1900);
        
        % amputation with 'always active' KCs:
        % number of useless always active KCs= 9(1-CL)*100
        
        KCs_alwaysActiv=round((9*(1-CL_rem_sparseKCs))*100);
        Homog_KCs_alwaysActiv=round((9*(1-CL_Homog_rem_sparseKCs))*100);
        KCs_silent=100-(KCs_alwaysActiv);
        Homog_KCs_silent=100-Homog_KCs_alwaysActiv;
        ValueMax=max(Y(:));
        ValueMax_H=max(YHomog_Ftheta(:));
        Y(i1(end-99:end),:)=[];
        YHomog_Ftheta(i2(end-99:end),:)=[];
        Y(end+1:end+KCs_alwaysActiv,:)=repmat(ValueMax.*ones(1,odors*numTrials),KCs_alwaysActiv,1);
        YHomog_Ftheta(end+1:end+Homog_KCs_alwaysActiv,:)=repmat(ValueMax_H.*ones(1,odors*numTrials),Homog_KCs_alwaysActiv,1);
        
        % amputation with 'silent' KCs: 
        Y(end+1:end+KCs_silent,:)=repmat(zeros(1,odors*numTrials),KCs_silent,1);
        YHomog_Ftheta(end+1:end+Homog_KCs_silent,:)=repmat(zeros(1,odors*numTrials),Homog_KCs_silent,1);
       
         
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
       
             % get the ang. distance after discarding the most specialised KCs            
             S_PW_AngDist(:,:,randomTrials)=AngDist_Pw_Calc(Y);
             H_PW_AngDist(:,:,randomTrials)=AngDist_Pw_Calc(YHomog_Ftheta);
       
end

% save the models performances & angular distances
save( strcat(['Excl_mostSpecificKCs_rand_homogModels_',num2str(odors),'odors.mat']),'test_p_ra','test_p_raH_FixedTheta','S_PW_AngDist','H_PW_AngDist');


