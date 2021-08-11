%% Date: 11 August 2021
% author: Nada Abdelrahman
%  this code replicates the results in Fig2E,3(A2,C-E),
%  and S2(B,D-G) of:
%  "Compensatory variability in network parameters
%  enhances memory performance in the Drosophila mushroom body".
%  Highlighting that inter-KCs variability is beneficial 
%  in non sparse (dense) coding levels regime. 
%%

clear all
clc
load('hallem_olsen.mat');
PN = hallem_olsen(1:110,:)'; 
numtrainingSamples=15; %% artificial plus linearly dependent responses 

% uncomment this to run the analysis for 20 odors task.
% mulOd=20;
mulOd=100;

lrs=10;
LRs=10.^[-5 -4 -3 -2.75 -2.5 -2.25 -2 -1 0 1]; 
Crange=1;           
odors=mulOd;
numTrials = 30;
n =2000; % number of neurons in the hidden layer
m=24; 
k=0;
%% create artificial odors, n odors
for Pn=1:24
    [prob,bins]=hist(PN(Pn,:),100);
    prob=prob/sum(prob);
    % uncomment the next line if you want to use new set of odors responses
    % and comment lines (37 to 81)
    %PNs_1(Pn,k+1:k+mulOd)=randsample(bins,mulOd,'true',prob);
    
    binsAll(Pn,:) = bins; 
    
end

%% these are the saved PN responses we used through our simulations:
% comment this section if you want to use new set of odors
load('Data_submitted_fly_wNoise11.mat','PNtrials');

% Pre-scaled PNs responses: the base 100 fictitious odors.
PNs_1_=PNtrials(:,:,1);

%% preprocessing step: recover the rescaling factors
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
PNs_1 = (PNs_1_ - p(2))/p(1);

%% these are the saved PN responses we used in our simulations:
% comment this section if you want to use new set of odors 

PNs=PNs_1;
% uncomment this to use 20 odors subset of the 100 odors:
% Subset_RandOds=randsample(100,20);
% PNs=PNs_1(:,Subset_RandOds);

C_SoftMax=1;
x=PNs;
lifeTime_Spar=zeros(n,50,2);
H_lifeTime_Spar=zeros(n,50,2);
S_PW_AngDist=zeros(odors,odors,50,2);
H_PW_AngDist=zeros(odors,odors,50,2);

for randomTrials=1:50

    l=1;
    l2=1;
    classAction1=randsample([1:odors],round(odors/2),'false'); 
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
        

    for k=1:n

        for j=1:length(PnToKc{k})

          whichPN = PnToKc{k}(j);

          % pick random weight from a log normal distribution that
            % roughtly fits the Turner distribution

           thisWeight = exp(-0.0507+0.3527*randn(1));


          thisW(whichPN, k) = thisW(whichPN, k) + thisWeight;

        end
    end


    for k=1:n

        for j=1:length(HomogPnToKc{k})


          whichPN_homog= HomogPnToKc{k}(j);

          thisWeightHomo=1; %% homogenous equal unity weights connecting KCs to PNs.

          thisW_HomogModel(whichPN_homog,k)= thisWeightHomo+ thisW_HomogModel(whichPN_homog,k); 


        end
    end


    noise=1;
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


     PNtrials(PNtrials<0)=0;
     PNtrials=rescale(PNtrials,0,5);
                        
    %% tuning starts: tune for CL=10% train and test the network, then for 30%, 70% & 90%, train and test this network. 
    %  Here APL=0
    CLConds=[0.1,0.3,0.5,0.7,0.9];
    CLConds_noAPL=[0.1,0.3,0.5,0.7,0.9];

    for CLCond=1:size(CLConds,2)

                tic
                APLgainP= zeros(1,2);  
                C_thetaS=1;
                C_thetaH=1;

                CL=CLConds(CLCond);
                CLnoAPL=CLConds_noAPL(CLCond);
%% 
                % Tuning the Random model:
                
                T=15;
                thetaS=abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
                thetaS(thetaS>70)=70;
                thetaS(thetaS<0.01)=0.01;
                constraints=0; 
                A=zeros(n,odors*numtrainingSamples);
                Y=zeros(n,odors*numtrainingSamples);
                eta=0.1; 

                while(~constraints)

                % with inhibition gain absent, scaling thetas distribution to acheive 
                % the sparsity constraint= 20% when APL feedback is inhibited.
                for trial = 1:(odors*numtrainingSamples) 

                    A(:,trial) = thisW'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)-(C_thetaS.*thetaS))>0 ).*( A(:,trial)-(C_thetaS.*thetaS));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);

                end
                InhAbs_mSimp=mean(codingLevelDummy);    
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(thetaS,1,odors*numtrainingSamples));
                Grad= ((InhAbs_mSimp)-CLnoAPL)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
                C_thetaS=C_thetaS - eta.*(Grad);

                if (C_thetaS<0)
                error('the scale factor in the random model is -ve')
                end
    
                % check if the constraints are satisfied
                for trial = 1:(odors*numtrainingSamples) 

                    Activations_S(:,trial) = thisW'*PNtrials(:,trial);
                    Y_S(:,trial)=(( Activations_S(:,trial)-(C_thetaS.*thetaS))>0 ).*( Activations_S(:,trial)-(C_thetaS.*thetaS));
                    codingLevelDummy(trial)=  (sum(Y_S(:,trial)>0,1)/n); 
                end

                 InhAbs_CL=mean(codingLevelDummy);

                for trial = 1:(odors*numtrainingSamples) 

                    Activation(:,trial) = thisW'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgainP(1))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS))>0 ).*( Activation(:,trial)-APLgainP(1)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaS.*thetaS));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n); 

                end
                 CL_=mean(codingLevelDummy);
                 
                 constraints= (abs(InhAbs_CL-CLnoAPL)<0.01) & (abs(CL_-CL)<0.01);

                end

                CLevelP(1)=CL_;
                thetaS=(C_thetaS.*thetaS);
                INHAbs_CLP(1)=InhAbs_CL;

%% 
                % Tuning the homogenous model:
                
                constraints=0; 
                thetaH_Ftheta=5+rand(1);
                A=zeros(n,odors*numtrainingSamples);
                Y=zeros(n,odors*numtrainingSamples);
                HomoFtheta_Inh=0;

                while(~constraints)

                % with inhibition gain absent, scaling thetas distribution to acheive 
                % the sparsity constraint= 20% when APL feedback is inhibited.    
                eta=0.1;   
                for trial = 1:(odors*numtrainingSamples) 

                    A(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y(:,trial)=(( A(:,trial)-(C_thetaH*thetaH_Ftheta))>0 ).*( A(:,trial)-(C_thetaH*thetaH_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n);

                end
                InhAbs_mHomog=mean(codingLevelDummy);    
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaH_Ftheta,n,odors*numtrainingSamples));
                Grad= ((InhAbs_mHomog)-CLnoAPL)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                C_thetaH=C_thetaH- (eta.*Grad);

                if (C_thetaH<0)
                   error('the scale factor in black model is -ve')
                end

                
                % check if the constraints are satisfied 
                for trial = 1:(odors*numtrainingSamples) 

                    Activations_H(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y_H(:,trial)=(( Activations_H(:,trial)-(C_thetaH*thetaH_Ftheta))>0 ).*( Activations_H(:,trial)-(C_thetaH*thetaH_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y_H(:,trial)>0,1)/n); 
                end

                 InhAbs_CL=mean(codingLevelDummy);

                for trial = 1:(odors*numtrainingSamples) 

                    Activation(:,trial) = thisW_HomogModel'*PNtrials(:,trial);
                    Y(:,trial)=(( Activation(:,trial)-(APLgainP(2))*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH*thetaH_Ftheta))>0 ).*( Activation(:,trial)-APLgainP(2)*repmat(sum(Activation(:,trial),1),n,1)-(C_thetaH*thetaH_Ftheta));
                    codingLevelDummy(trial)=  (sum(Y(:,trial)>0,1)/n); 

                end
                 CL_=mean(codingLevelDummy);

                constraints= (abs(InhAbs_CL-CLnoAPL)<0.01) & (abs(CL_-CL)<0.01);

                end

               CLevelP(2)=CL_;
               INHAbs_CLP(2)=InhAbs_CL;
               thetaH_Ftheta=(C_thetaH*thetaH_Ftheta);

               toc
%% 
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

              if (CL==0.1 || CL==0.9)                
                  %% LT sparsity
                   Kpats=odors*numTrials;

                   Spar=(1/(1-(1/Kpats)))*(1-(( sum( (Y./Kpats),2 ) ).^2./(sum( ( (Y.^2)./Kpats),2 ) )));

                   HomogenousSpar= (1/(1-(1/Kpats)))*(1-((sum( (YHomog_Ftheta./Kpats),2)).^2./(sum( ( (YHomog_Ftheta.^2)./Kpats),2 ))));

                   lifeTime_Spar(:,randomTrials,l)= Spar;

                   H_lifeTime_Spar(:,randomTrials,l)= HomogenousSpar;

                    % standard deviation in the lifetime spar.values 
                    H_LTSpar(randomTrials,l)=std(HomogenousSpar(~isnan(HomogenousSpar)));   
                    S_LTSpar(randomTrials,l)=std(Spar(~isnan(Spar)));

                    % angular distance PW 
                     S_PW_AngDist(:,:,randomTrials,l)=AngDist_Pw_Calc(Y);
                     H_PW_AngDist(:,:,randomTrials,l)=AngDist_Pw_Calc(YHomog_Ftheta);

                    l=l+1;

              end


              YHomog_Fthetatemp1=reshape(YHomog_Ftheta,n,odors,numTrials);
              Ytemp1= reshape(Y,n,odors,numTrials);
              
              % rescale the KC responses to be from 0-1
              Ytemp=rescale(Ytemp1);
              YHomog_Fthetatemp=rescale(YHomog_Fthetatemp1);
              YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
              Ytr=Ytemp(:,:,1:numtrainingSamples);

              % calculate DBI between odor pairs:
              
              HF_DBi_odoursPairs(:,:,randomTrials,CLCond)= dbiCalc(YHomog_Fthetatemp);
              S_DBi_odoursPairs (:,:,randomTrials,CLCond)= dbiCalc(Ytemp);

            for l_r=1:lrs  

                WopAllOdours=1*rand(n,2);              
                WopAllOdoursHomog_Ftheta=WopAllOdours;
                alpha=LRs(l_r);
                c=1;
                ceq=1;
                ch=1;
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

                 end

                 %% Testing the tuned models:
                 for c=1:Crange
                     C=C_SoftMax*(10^c);

                     [acc]=testingRandomModel_accuracy(C,WopAllOdours,classAction1,numTrials,numtrainingSamples,Ytemp);

                    test_p_ra(randomTrials, l_r,CLCond)=acc;


                  [accH1]=HomogenousModel_KernelTesting(C,WopAllOdoursHomog_Ftheta,PNtrials, thetaH_Ftheta,classAction1,numTrials,numtrainingSamples,YHomog_Fthetatemp);

                    test_p_raH_FixedTheta(randomTrials, l_r,CLCond)=accH1;  

                 end

            end
            %% calculate angular distance between 'bad' and 'good'
               %odors classes: Tensor [ random fly networks (20) x 1]

              odSet= [1:odors];
              class2= odSet(~ismember(odSet,classAction1));                   
              labels_acrossTrials= repmat( odSet',numTrials,1);

               S1Resp= Ytemp(:,classAction1,:);
               S1RespAcTrials=reshape(S1Resp,n,(odors/2)*numTrials);
               S1RespAcTrialsCentroid=mean(S1RespAcTrials,2);
               S1RespVar=sqrt(mean(vecnorm((S1RespAcTrials-S1RespAcTrialsCentroid ) )));

               S2Resp= Ytemp(:,class2,:);
               S2RespAcTrials=reshape(S2Resp,n,(odors/2)*numTrials);
               S2RespAcTrialsCentroid=mean(S2RespAcTrials,2);
               S2RespVar=sqrt(mean(vecnorm((S2RespAcTrials-S2RespAcTrialsCentroid) )));

               S1RespHo_Ftheta= YHomog_Fthetatemp(:,classAction1,:);
               S1RespAcTrialsHo_Ftheta=reshape(S1RespHo_Ftheta,n,(odors/2)*numTrials);
               S1RespAcTrialsHoCentroid_Ftheta=mean(S1RespAcTrialsHo_Ftheta,2);
               S1RespVarHo_Ftheta= sqrt(mean(vecnorm(S1RespAcTrialsHo_Ftheta-S1RespAcTrialsHoCentroid_Ftheta)));

               S2RespHo_Ftheta= YHomog_Fthetatemp(:,class2,:);
               S2RespAcTrialsHo_Ftheta=reshape(S2RespHo_Ftheta,n,(odors/2)*numTrials);
               S2RespAcTrialsHoCentroid_Ftheta=mean(S2RespAcTrialsHo_Ftheta,2);
               S2RespVarHo_Ftheta= sqrt(mean(vecnorm((  S2RespAcTrialsHo_Ftheta- S2RespAcTrialsHoCentroid_Ftheta) )));  

             % calculate the angular distance between 'good' and 'bad'
             % odors clusters:
             S_Cw_AngDist(randomTrials,CLCond)= (acos(( (S1RespAcTrialsCentroid'/norm(S1RespAcTrialsCentroid'))*(S2RespAcTrialsCentroid/norm(S2RespAcTrialsCentroid))) )/ (0.5*pi) ) ;
             HF_Cw_AngDist(randomTrials,CLCond)= (acos(( (S1RespAcTrialsHoCentroid_Ftheta'/norm(S1RespAcTrialsHoCentroid_Ftheta'))*(S2RespAcTrialsHoCentroid_Ftheta/norm(S2RespAcTrialsHoCentroid_Ftheta))) )/ (0.5*pi) ) ;
             
             % calculate the DBI between 'good' and 'bad' odors clusters:
             S_Cw_DBI(randomTrials,CLCond)=  (S1RespVar+S2RespVar)/( norm( S1RespAcTrialsCentroid-S2RespAcTrialsCentroid));
             HF_Cw_DBI(randomTrials,CLCond)= (S1RespVarHo_Ftheta+S2RespVarHo_Ftheta)/(norm(S1RespAcTrialsHoCentroid_Ftheta- S2RespAcTrialsHoCentroid_Ftheta));

             if (CL==0.1 || CL==0.9)
                 
                 % calculate KCs' valence specificity: 
                 V_spec_R(:,randomTrials,l2)= abs( sum(Ytemp1(:,classAction1,:),[2 3])- sum(Ytemp1(:,class2,:), [2 3]) )./ sum(Ytemp1,[2 3]); 
                 V_spec_Homog(:,randomTrials,l2)= abs( sum(YHomog_Fthetatemp1(:,classAction1,:),[2 3])- sum(YHomog_Fthetatemp1(:,class2,:), [2 3]) )./ sum(YHomog_Fthetatemp1,[2 3]); 
                 
                 l2=l2+1;
             end
       % save the tuned fly models to the current directory:
       save(strcat([pwd,'CL_performance_random_homog',num2str(odors),'_',num2str(CL),'_',num2str(randomTrials),'.mat']),'PNtrials','classAction1','thisW','thisW_HomogModel','APLgainP', 'thetaH_Ftheta','thetaS');

    end              
       
end
% save the angular distance, lifetime sparsity, and Valence specificity
% measures for both fly models, under the different coding levels: CL=0.1 & 0.9 (without APL)
save(strcat('CL_APL_Iszero_performance_separabilityMetrics_',num2str(odors),'.mat'),'test_p_ra','test_p_raH_FixedTheta','S_Cw_AngDist','HF_Cw_AngDist','S_PW_AngDist','H_PW_AngDist','lifeTime_Spar','S_LTSpar','H_lifeTime_Spar','H_LTSpar'...
            ,'HF_DBi_odoursPairs','S_DBi_odoursPairs','HF_Cw_DBI','S_Cw_DBI','V_spec_R','V_spec_Homog');
