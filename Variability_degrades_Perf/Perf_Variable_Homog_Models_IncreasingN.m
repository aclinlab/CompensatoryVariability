%% Date:11 August 2021
% Author: Nada Abdelrahman
% this code replicates the results of Fig2D
% in: "Compensatory variability in network parameters
% enhances memory performance in the Drosophila mushroom body"
% Inter-KC variability increases odors encoding separability
% as the number of PN inputs/KC (N) increases
%%

clear all
clc
load('hallem_olsen.mat');
numtrainingSamples=15; %% artificial plus linearly dependent responses 
mulOd=100;
lrs=10;
Crange=1;           
odors=mulOd;
numTrials = 30;
n =2000; % number of neurons in the hidden layer
m=24; 
LRs=10.^[-5 -4 -3 -2.75 -2.5 -2.25 -2 -1 0 1]; 
C_SoftMax=1;
NScales=1;

% uncomment any of the following N values 
% you wish to carry the analysis for:

NPNs=24;
% NPNs=18;
% NPNs=12;
% NPNs=6;

for randomTrials=1:30

    % load the same odor inputs and valence assignments in Fig2B&C
    load(strcat(['2B_C_VarDegradesPerf_fly_wNoise',num2str(randomTrials),num2str(2)]),'PNtrials','classAction1');
    
    clawsNo=ones([1,n])*NPNs; 
    HomogClaws= ones([1,n])*NPNs; 
    HomogPNsperKC= HomogClaws;
       
    for i=1:n
        HomogPnToKc{i}= randsample(m, HomogPNsperKC(i));
    end % random initilaization of the weights

    %initialize the weights matrix between KC and PN

    thisW = zeros(m, n);
    thisW_HomogModel=zeros(m,n);

    for k=1:n

        for j=1:length(HomogPnToKc{k})

          whichPN = HomogPnToKc{k}(j);

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

%%    tune the random fly model for the sparsity constraints:
                
            tic
            APLgainP= zeros(1,3);  
            C_thetaS=1;
            C_thetaH=1;
            CL=0.1;
            CLnoAPL=0.2;

            T=15;
            thetaS=abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
            thetaS(thetaS<0.01)=0.01;
            constraints=0; 
            eta=0.1; 

            while(~constraints)

            A=multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
            Y=(( A-(C_thetaS.*thetaS))>0 ).*( A-(C_thetaS.*thetaS));
            Y=reshape(Y,n,odors*numtrainingSamples);
            codingLevelDummy= sum(Y>0)./n;

            InhAbs_mSimp=mean(codingLevelDummy);    
            depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
            depsi1_dy(isnan(depsi1_dy))=0;
            depsi1_dtheta= -(Y>0).* depsi1_dy.* (repmat(thetaS,1,odors*numtrainingSamples));
            Grad= ((InhAbs_mSimp)-CLnoAPL)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:)));
            C_thetaS=C_thetaS - eta.*(Grad);

            if (C_thetaS<0)
            error('the scale factor in the random model is -ve')
            end

%           replicating the sparsity level of the KCs in real MB network
%           CL=10%  
            eta_2=0.00000001;

            ActivationsEqualizeddummy = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
            YEqualizeddummy=(( ActivationsEqualizeddummy-(APLgainP(1))*repmat(sum(ActivationsEqualizeddummy,1),n,1)-(C_thetaS.*thetaS))>0 ).*( ActivationsEqualizeddummy-APLgainP(1)*repmat(sum(ActivationsEqualizeddummy,1),n,1)-(C_thetaS.*thetaS));
            YEqualizeddummy=reshape(YEqualizeddummy,n,odors*numtrainingSamples);
            codingLevelEqualizedDummy=  (sum(YEqualizeddummy>0)./n);

            CLRand=mean(codingLevelEqualizedDummy);
            dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
            dsig_dy(isnan(dsig_dy))=0;
            dAct_dalpha= sum(  reshape(ActivationsEqualizeddummy,n,odors*numtrainingSamples) ,1);
            dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
            Grad_alpha= ((CLRand)-CL)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
            APLgainP(1)= APLgainP(1)- eta_2*(Grad_alpha);

            if (APLgainP(1)<0)
               error('APL gain 1 <0')
            end

            % check if the constraints are satisfied
            Activations_S = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
            Y_S=(( Activations_S-(C_thetaS.*thetaS))>0 ).*( Activations_S-(C_thetaS.*thetaS));
            Y_S=reshape(Y_S,n,odors*numtrainingSamples);
            codingLevelDummy=  (sum(Y_S>0)./n); 
            
            InhAbs_CL=mean(codingLevelDummy);

            Activation = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
            Y=(( Activation-(APLgainP(1))*repmat(sum(Activation,1),n,1)-(C_thetaS.*thetaS))>0 ).*( Activation-APLgainP(1)*repmat(sum(Activation,1),n,1)-(C_thetaS.*thetaS));
            Y=reshape(Y,n,odors*numtrainingSamples);
            codingLevelDummy=  (sum(Y>0)./n); 

            CL_=mean(codingLevelDummy);


            constraints= (abs(InhAbs_CL-CLnoAPL)<0.01) & (abs(CL_-CL)<0.01);

            end

            CLevelP(1)=CL_;
            thetaS=(C_thetaS.*thetaS);
            INHAbs_CLP(1)=InhAbs_CL;


%%     tune the homogeneous model: 
                
            constraints=0; 
            thetaH_Ftheta=20+rand(1);             
            HomoFtheta_Inh=0;

            while(~constraints)

            % with inhibition gain absent, scaling thetas distribution to acheive 
            % the sparsity constraint= 20% when APL feedback is inhibited.    
            eta=0.1;   
            A = multiprod(thisW_HomogModel',PNtrials(:,:,1:numtrainingSamples));
            Y=(( A-(C_thetaH*thetaH_Ftheta))>0 ).*( A-(C_thetaH*thetaH_Ftheta));
            Y=reshape(Y,n,odors*numtrainingSamples);
            codingLevelDummy=  (sum(Y>0)./n);

            InhAbs_mHomog=mean(codingLevelDummy);    
            depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
            depsi1_dy(isnan(depsi1_dy))=0;
            depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaH_Ftheta,n,odors*numtrainingSamples));
            Grad= ((InhAbs_mHomog)-CLnoAPL)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
            C_thetaH=C_thetaH- (eta.*Grad);

            if (C_thetaH<0)
               error('the scale factor in black model is -ve')
            end

            % replicating the sparsity level of the KCs in real MB network
            % CL=10%        
            eta_2=0.000000001;
            ActivationsEqualizeddummy = multiprod(thisW_HomogModel',PNtrials(:,:,1:numtrainingSamples));
            YEqualizeddummy=(( ActivationsEqualizeddummy-(APLgainP(2))*repmat(sum(ActivationsEqualizeddummy,1),n,1)-(C_thetaH*thetaH_Ftheta))>0 ).*( ActivationsEqualizeddummy-APLgainP(2)*repmat(sum(ActivationsEqualizeddummy,1),n,1)-(C_thetaH*thetaH_Ftheta));
            YEqualizeddummy=reshape(YEqualizeddummy,n,odors*numtrainingSamples);
            codingLevelEqualizedDummy=  (sum(YEqualizeddummy>0)./n);

            CLAllFixed=mean(codingLevelEqualizedDummy);
            dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
            dsig_dy(isnan(dsig_dy))=0;
            dAct_dalpha= sum( reshape(ActivationsEqualizeddummy,n,odors*numtrainingSamples) ,1);
            dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
            Grad_alpha= ((CLAllFixed)-CL)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
            APLgainP(2)= APLgainP(2)- eta_2*(Grad_alpha);

            if (APLgainP(2)<0)
               error('APL gain 2 <0')
            end
            % check if the constraints are satisfied 
            Activations_H = multiprod(thisW_HomogModel',PNtrials(:,:,1:numtrainingSamples));
            Y_H=(( Activations_H-(C_thetaH*thetaH_Ftheta))>0 ).*( Activations_H-(C_thetaH*thetaH_Ftheta));
            Y_H=reshape(Y_H,n,odors*numtrainingSamples);
            codingLevelDummy=  (sum(Y_H>0)./n); 
            
            InhAbs_CL=mean(codingLevelDummy);


            Activation = multiprod(thisW_HomogModel',PNtrials(:,:,1:numtrainingSamples));
            Y=(( Activation-(APLgainP(2))*repmat(sum(Activation,1),n,1)-(C_thetaH*thetaH_Ftheta))>0 ).*( Activation-APLgainP(2)*repmat(sum(Activation,1),n,1)-(C_thetaH*thetaH_Ftheta));
            Y=reshape(Y,n,odors*numtrainingSamples);
            codingLevelDummy=  (sum(Y>0)./n); 

            CL_=mean(codingLevelDummy);

            constraints= (abs(InhAbs_CL-CLnoAPL)<0.01) & (abs(CL_-CL)<0.01);
            end

           CLevelP(2)=CL_;
           INHAbs_CLP(2)=InhAbs_CL;
           thetaH_Ftheta=(C_thetaH*thetaH_Ftheta);

%%         tune the varying w, fixed N, fixed theta model:
           
         thetaS_Ftheta= 30+rand(1);
         constraints=0; 
                      
        C_thetaS_Ftheta=1;
        while(~constraints)

                eta=0.1;   
                A = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
                Y=(( A- (C_thetaS_Ftheta*thetaS_Ftheta) )>0 ).*( A-(C_thetaS_Ftheta*thetaS_Ftheta));
                Y=reshape(Y,n,odors*numtrainingSamples);
                codingLevelDummy=  (sum(Y>0)./n);
                InhAbs_mSimp=mean(codingLevelDummy);    
        
                %% we want to change the theta so to achieve InhAbs_CL=2xCL
                depsi1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                depsi1_dy(isnan(depsi1_dy))=0;
                depsi1_dtheta= -(Y>0).* depsi1_dy.*(repmat(thetaS_Ftheta,n,odors*numtrainingSamples));
                Grad= ((InhAbs_mSimp)-0.20)*(1/(n*odors*numtrainingSamples)).*(sum(depsi1_dtheta(:) ));
                C_thetaS_Ftheta=C_thetaS_Ftheta - eta.*(Grad);

                   %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%            
                eta_2=0.00000001;
                ActivationsEqualizeddummy = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
                YEqualizeddummy=(( ActivationsEqualizeddummy-(APLgainP(3))*repmat(sum(ActivationsEqualizeddummy,1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta))>0 ).*( ActivationsEqualizeddummy-APLgainP(3)*repmat(sum(ActivationsEqualizeddummy,1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta));
                YEqualizeddummy=reshape(YEqualizeddummy,n,odors*numtrainingSamples);
                codingLevelEqualizedDummy=  (sum(YEqualizeddummy>0)./n);
                CLRandFixedTh=mean(codingLevelEqualizedDummy);
                dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
                dsig_dy(isnan(dsig_dy))=0;
                dAct_dalpha= sum(reshape(ActivationsEqualizeddummy,n,odors*numtrainingSamples),1);
                dsig_dalpha= -(YEqualizeddummy>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                Grad_alpha= ((CLRandFixedTh)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
                APLgainP(3)= APLgainP(3)- eta_2*(Grad_alpha);

            %% check if the constraints are satisfied 
                Activations_SF = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
                Y_SF=(( Activations_SF-(C_thetaS_Ftheta*thetaS_Ftheta))>0 ).*( Activations_SF-(C_thetaS_Ftheta*thetaS_Ftheta));
                Y_SF=reshape(Y_SF,n,odors*numtrainingSamples);
                codingLevelDummy=  (sum(Y_SF>0)./n); 
             InhAbs_CL=mean(codingLevelDummy);

                Activation = multiprod(thisW',PNtrials(:,:,1:numtrainingSamples));
                Y=(( Activation-(APLgainP(3))*repmat(sum(Activation,1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta))>0 ).*( Activation-APLgainP(3)*repmat(sum(Activation,1),n,1)-(C_thetaS_Ftheta*thetaS_Ftheta));
                Y=reshape(Y,n,odors*numtrainingSamples);
                codingLevelDummy=  (sum(Y>0)./n); 
             CL_=mean(codingLevelDummy);

             constraints= (abs(InhAbs_CL-CLnoAPL)<0.01) & (abs(CL_-CL)<0.01);

        end
        thetaS_Ftheta=(C_thetaS_Ftheta*thetaS_Ftheta);
      
        %%
        toc

       Activations=zeros(n,odors*numTrials);
       ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numTrials);
       Activations_Ftheta=zeros(n,odors*numTrials);

       Y=zeros(n,odors*numTrials);
       YHomog_Ftheta=zeros(n,odors*numTrials);
       Y_Ftheta=zeros(n,odors*numTrials);

       for trial = 1:(odors*numTrials)


           ActivationsHomogenousdummy_Ftheta(:,trial) = thisW_HomogModel'*PNtrials(:,trial  );
           YHomog_Ftheta(:,trial)=(( ActivationsHomogenousdummy_Ftheta(:,trial)-(APLgainP(2))*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomogenousdummy_Ftheta(:,trial)-APLgainP(2)*repmat(sum(ActivationsHomogenousdummy_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);                   

           Activations(:,trial) = thisW'*PNtrials(:,trial );
           Y(:,trial)=(( Activations(:,trial)-(APLgainP(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgainP(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);

           Activations_Ftheta(:,trial) = thisW'*PNtrials(:,trial );
           Y_Ftheta(:,trial)=(( Activations_Ftheta(:,trial)-(APLgainP(3))*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta)>0 ).*( Activations_Ftheta(:,trial)-APLgainP(3)*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta);


       end

        for l_r=1:lrs  
               
            WopAllOdours=1*rand(n,2);              
            WopAllOdoursHomog_Ftheta=WopAllOdours;
            WopAllOdours_Ftheta=WopAllOdours;
            alpha=LRs(l_r);                
            c=1;
            ceq=1;
            ch=1;
            YHomog_Fthetatemp=reshape(YHomog_Ftheta,n,odors,numTrials);
            Y_Fthetatemp= reshape(Y_Ftheta,n,odors,numTrials); 
            Ytemp= reshape(Y,n,odors,numTrials);
            
            % rescale the KC responses to be from 0-1
             Ytemp=rescale(Ytemp);
             YHomog_Fthetatemp=rescale(YHomog_Fthetatemp);
             Y_Fthetatemp=rescale(Y_Fthetatemp);
             YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
             Y_Fthetatr=Y_Fthetatemp(:,:,1:numtrainingSamples);  
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


                    if( ~ isempty(find(classAction1==odour)) )              
                         delta =  exp( -(alpha/mean(Y_Fthetatr(:)))* sum(Y_Fthetatr(:,odour,:),3) );
                         WopAllOdours_Ftheta(:,2)= WopAllOdours_Ftheta(:,2) .*delta;

                    else
                          delta = exp(- (alpha/mean(Y_Fthetatr(:)))* sum(Y_Fthetatr(:,odour,:),3) );
                          WopAllOdours_Ftheta(:,1)= WopAllOdours_Ftheta(:,1) .*delta;
                    end


             end
                
             % testing the tuned models:
                 
             for c=1:Crange
                 
                
                 C=C_SoftMax*(10^c);
               
                 [acc]=testingRandomModel_accuracy(C,WopAllOdours,classAction1,numTrials,numtrainingSamples,Ytemp);
                
                 test_p_ra(randomTrials, l_r)=acc;
                
                 [acc2]=testingVaryWeightsModel_accuracy(C,WopAllOdours_Ftheta,classAction1,numTrials,numtrainingSamples,Y_Fthetatemp);
                
                 test_p_ra_Ftheta(randomTrials, l_r)=acc2;
                
                
                 [accH1]=HomogenousModel_KernelTesting(C,WopAllOdoursHomog_Ftheta,PNtrials, thetaH_Ftheta,classAction1,numTrials,numtrainingSamples,YHomog_Fthetatemp);

                 test_p_raH_FixedTheta(randomTrials, l_r)=accH1;  
                
 
      
             end
      
        end
        
end
% save models accuracies for this condition of N PN inputs:
save(strcat(['2D_Perf_',num2str(NPNs),'N.mat']),'test_p_ra','test_p_ra_Ftheta','test_p_raH_FixedTheta');