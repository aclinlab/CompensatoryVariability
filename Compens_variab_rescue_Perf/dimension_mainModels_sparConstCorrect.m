%% code for dimension calculation for all the tuned models
%% the main modeels, coding level constraint corrected by CTheta
clear all
clc
CurrDir=pwd;
numtrainingSamples=15;
numTrials = 30;
mulOd=1500;
load('hallem_olsen.mat');

PN = hallem_olsen(1:110,:)';
PNs=zeros(24,mulOd);
% get 1000 odors from the 110 PNs
k=0; 
        %% create artificial odors, n odors

            for Pn=1:24
                [prob,bins]=hist(PN(Pn,:),100);
                prob=prob/sum(prob);

                PNs(Pn,k+1:k+mulOd)=randsample(bins,mulOd,'true',prob);

            end
     noise=1*.1;
     PNfulltrials = zeros(24, mulOd, numTrials);

     PNfulltrials(:,:,1) =PNs;

        for t = 1:numTrials-1
                            PNfulltrials(:,:,t+1) = PNs + ...
                                getPNStdevBhandawat(PNs) .* ...
                                noise.*randn(24, mulOd);
        end

PNfulltrials(PNfulltrials<0)=0;
PNfulltrials=rescale(PNfulltrials,0,5);
n=2000;
ods=10;
ll=20;
% ll=50;
CL=1;
% for CL=1:1
%     CLs=[0.9];
    
for fly=1:ll

 load(strcat(CurrDir,'/TunedFlies_allModels_oldOdors_with_round',['_fly_wNoise',num2str(fly),num2str(1),'.mat']) ); 
 odors= mulOd*(ods/10);
 x=PNfulltrials(:,1:odors,:);
 PNtrials=reshape(x,24,odors*numTrials);
         
   Activations=zeros(n,odors*numTrials);
   ActivationsEqualized=zeros(n,odors*numTrials);
   ActivationsHomogenousdummy_Ftheta=zeros(n,odors*numTrials);
   Activations_comp2= zeros(n,odors*numTrials);
   Activations_Kenn= zeros(n,odors*numTrials);
   Activations_theta_activity_hoemeo= zeros(n,odors*numTrials);
   Activations_inhibPlast=zeros(n,odors*numTrials);
   Activations_comp2_noxjk=zeros(n,odors*numTrials);
   Activations_comp2_wHy= zeros(n,odors*numTrials);

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
                    Y_theta_activity_homeo(:,trial)=(( Activations_theta_activity_hoemeo(:,trial)-(APLgains(5))*repmat(sum(Activations_theta_activity_hoemeo(:,trial),1),n,1)-theta_Activity_homeo)>0 ).*( Activations_theta_activity_hoemeo(:,trial)-APLgains(5)*repmat(sum(Activations_theta_activity_hoemeo(:,trial),1),n,1)-theta_Activity_homeo);
        
       
                   
                   Activations(:,trial) = thisW'*PNtrials(:,trial );
                    Y(:,trial)=(( Activations(:,trial)-(APLgains(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgains(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
       
                   

                    ActivationsEqualized(:,trial) = thisW_equalizedModel'*PNtrials(:,trial );
                    YEqualized(:,trial)=(( ActivationsEqualized(:,trial)-(InhibitionGain)*repmat(sum(ActivationsEqualized(:,trial),1),n,1)-theta)>0 ).*( ActivationsEqualized(:,trial)-InhibitionGain*repmat(sum(ActivationsEqualized(:,trial),1),n,1)-theta);
        
        
                    Activations_comp2(:,trial) = thisW_ActivityBasedComp'*PNtrials(:,trial );
                    Y_comp2(:,trial)=(( Activations_comp2(:,trial)-(APLgains(3) )*repmat(sum(Activations_comp2(:,trial),1),n,1)-theta_comp2)>0 ).*( Activations_comp2(:,trial)-APLgains(3)*repmat(sum(Activations_comp2(:,trial),1),n,1)-theta_comp2);
        
                    
                    Activations_Kenn(:,trial) = thisW'*PNtrials(:,trial );
                    Y_Kenn(:,trial)=(( Activations_Kenn(:,trial)-(APLgains(4) )*repmat(sum(Activations_Kenn(:,trial),1),n,1)-theta_Kenn)>0 ).*( Activations_Kenn(:,trial)-APLgains(4)*repmat(sum(Activations_Kenn(:,trial),1),n,1)-theta_Kenn);
        
                    Activations_inhibPlast(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials(:,trial );
                    Y_inhibPlast(:,trial)=(( Activations_inhibPlast(:,trial)-(APLgains_model6 ).*repmat(sum(Activations_inhibPlast(:,trial),1),n,1)-theta_inhibitionPlast)>0 ).*( Activations_inhibPlast(:,trial)-APLgains_model6.*repmat(sum(Activations_inhibPlast(:,trial),1),n,1)-theta_inhibitionPlast);
                        
                   Activations_comp2_noxjk(:,trial) = thisW_ActivityBasedComp_noxjk'*PNtrials(:,trial );
                   Y_comp2_noxjk(:,trial)=(( Activations_comp2_noxjk(:,trial)-(APLgains_noxjk )*repmat(sum(Activations_comp2_noxjk(:,trial),1),n,1)-theta_comp2_noxjk)>0 ).*( Activations_comp2_noxjk(:,trial)-APLgains_noxjk*repmat(sum(Activations_comp2_noxjk(:,trial),1),n,1)-theta_comp2_noxjk);

                    Activations_comp2_wHy(:,trial) = thisW_ActivityBasedComp_wHy'*PNtrials(:,trial );
                    Y_comp2_wHy(:,trial)=(( Activations_comp2_wHy(:,trial)-(APLgains_wHy )*repmat(sum(Activations_comp2_wHy(:,trial),1),n,1)-theta_comp2_wHy)>0 ).*( Activations_comp2_wHy(:,trial)-APLgains_wHy*repmat(sum(Activations_comp2_wHy(:,trial),1),n,1)-theta_comp2_wHy);
 
                 
               end
               
  
               
               dim_S(fly,1)= dimInputCurrent(Calc_C(Y));
                
                dim_C(fly,1)=dimInputCurrent(Calc_C(YEqualized));
                
                dim_HF(fly,1)=dimInputCurrent(Calc_C(YHomogdummy_Ftheta));
                

                dim_C2(fly,1)=dimInputCurrent(Calc_C(Y_comp2));
                
                dim_C2_noxjk(fly,1)=dimInputCurrent(Calc_C(Y_comp2_noxjk));
%                 
                dim_C2_wHY(fly,1)=dimInputCurrent(Calc_C(Y_comp2_wHy));

                
                 dim_InhPlast(fly,1)=dimInputCurrent(Calc_C(Y_inhibPlast));
%                
                 dim_Kenn(fly,1)=dimInputCurrent(Calc_C(Y_Kenn));
                 dim_theta_Activity_homeo(fly,1)= dimInputCurrent(Calc_C(Y_theta_activity_homeo));
% %                  
                  
end
% end