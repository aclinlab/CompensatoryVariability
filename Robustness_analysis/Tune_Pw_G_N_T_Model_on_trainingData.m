function [thisW_equalizedModel_tune_is_train,theta_tune_is_train, InhibitionGain_tune_is_train,CL_, InhAbs_CL] = Tune_Pw_G_N_T_Model_on_trainingData(N_per_KC,PW_given_theta_and_n,W, thisW_equalizedModel_tune_is_train,...
                             PNtrials_tune_train, theta_tune_is_train,InhibitionGain_tune_is_train,odors,numtrainingSamples,C_theta_1)

% run the tuning but on the training set
constraints=0;
n=2000;                   
InhAbsTarget=0.20;
thisW_equalizedModel_0=thisW_equalizedModel_tune_is_train;
m=24;
    while(~constraints)
               
               eta=1;  
               
               for trial = 1:(odors*numtrainingSamples) 
                   
                    A_tune_tr(:,trial) = thisW_equalizedModel_tune_is_train'*PNtrials_tune_train(:,trial);
                    
                    Y_tune_tr(:,trial)=(( A_tune_tr(:,trial)- (C_theta_1.*theta_tune_is_train) )>0 ).*( A_tune_tr(:,trial)-(C_theta_1.*theta_tune_is_train));
                    codingLevel_eq(trial)=  (sum(Y_tune_tr(:,trial)>0,1)/n);

               end
               InhAbs_mComp_tune_tr=mean(codingLevel_eq);    
               %% we want to change the theta so to achieve InhAbs_CL=2xCL
               
              
               derr1_dy=(exp(0.9.*Y_tune_tr)./((1+exp(0.9.*Y_tune_tr)).^2));
               derr1_dy(isnan(derr1_dy))=0;
                  
               derr1_dtheta= -(Y_tune_tr>0).* derr1_dy.* (repmat(theta_tune_is_train,1,odors*numtrainingSamples));
              

               Grad= ((InhAbs_mComp_tune_tr)-InhAbsTarget)*(1/(n*odors*numtrainingSamples)).*(sum(derr1_dtheta(:)));
                   
               C_theta_1=C_theta_1 - eta.*(Grad);
               
               
              if (C_theta_1<0)
                     error('the scale factor in cyan model tune is train is -ve')
              end
              
               %% resample the weights
               thisW_equalizedModel_tune_is_train=zeros(m,n);

                for k=1:n

                 for j=1:m

                    if(thisW_equalizedModel_0(j,k))
                      %% sample the weights from the new fitted weights in the other script (modelling KC_PnWeights.m)

                      ThetaInd= round(((theta_tune_is_train(k)-0.01)/0.1)+1);
                      %% capping theta at 0.1..
                       ThetaInd(ThetaInd==0)=1;

                       this_KCWeights= PW_given_theta_and_n(N_per_KC(k)-1,ThetaInd,:);


                       thisWeight_equalizedModel= randsample(W,1,'true', this_KCWeights);

                      thisW_equalizedModel_tune_is_train(j,k)= thisW_equalizedModel_tune_is_train(j,k)+thisWeight_equalizedModel;
                    end


                 end
                end

                
               %% calculate the CL when there is inhibition
               eta_2=0.000001;
               
               for trial = 1:(odors*numtrainingSamples) 
                   
                    Activations_tutr(:,trial) = thisW_equalizedModel_tune_is_train'*PNtrials_tune_train(:,trial);
                    YEqualizedtutr(:,trial)=(( Activations_tutr(:,trial)-(InhibitionGain_tune_is_train)*repmat(sum(Activations_tutr(:,trial),1),n,1)-(C_theta_1.*theta_tune_is_train) )>0 ).*( Activations_tutr(:,trial)-InhibitionGain_tune_is_train*repmat(sum(Activations_tutr(:,trial),1),n,1)-(C_theta_1.*theta_tune_is_train));
                    codingLevel_eq_(trial)=  (sum(YEqualizedtutr(:,trial)>0,1)/n);

               end
               mComp=mean(codingLevel_eq_);
               
               dsig_dy_tutr=(exp(0.9.*YEqualizedtutr)./((1+exp(0.9.*YEqualizedtutr)).^2));
               dsig_dy_tutr(isnan(dsig_dy_tutr))=0;
              
               dAct_dalpha_= sum(Activations_tutr,1);
                   
               derr_dalpha= -(YEqualizedtutr>0).*(repmat(dAct_dalpha_,n,1)).*  dsig_dy_tutr;
                  

               Grad_alpha= ((mComp)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(derr_dalpha(:)));
               
               InhibitionGain_tune_is_train= InhibitionGain_tune_is_train- eta_2*(Grad_alpha);
                 
               if (InhibitionGain_tune_is_train<0)
                   
                   wh=1;
               end
                 
               
               %% check if the constraints are satisfied 
               Activations= zeros(n,odors*numtrainingSamples);
               Y= zeros(n,odors*numtrainingSamples);
               codingLevel_eq=zeros(1,odors*numtrainingSamples);
               
                for trial = 1:(odors*numtrainingSamples) 
                   
                    Activations(:,trial) = thisW_equalizedModel_tune_is_train'*PNtrials_tune_train(:,trial);
                    
                    Y(:,trial)=(( Activations(:,trial)-(C_theta_1.* theta_tune_is_train) )>0 ).*( Activations(:,trial)-(C_theta_1.* theta_tune_is_train));
                    
                    codingLevel_eq(trial)=  (sum(Y(:,trial)>0,1)/n); 
                end
                 
                 InhAbs_CL=mean(codingLevel_eq);
                 
               % 
               
               Activations= zeros(n,odors*numtrainingSamples);
               Y= zeros(n,odors*numtrainingSamples);
               codingLevel_eq=zeros(1,odors*numtrainingSamples);


                  for trial =1:(odors*numtrainingSamples) 
                   
                    Activations(:,trial) = thisW_equalizedModel_tune_is_train'*PNtrials_tune_train(:,trial);
                    Y(:,trial)=(( Activations(:,trial)-(InhibitionGain_tune_is_train)*repmat(sum(Activations(:,trial),1),n,1)-(C_theta_1.* theta_tune_is_train))>0 ).*( Activations(:,trial)-InhibitionGain_tune_is_train*repmat(sum(Activations(:,trial),1),n,1)-(C_theta_1.* theta_tune_is_train));
                    codingLevel_eq(trial)=  (sum(Y(:,trial)>0,1)/n); 

                  end
                 CL_=mean(codingLevel_eq);
                 
                
                 constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.1) &( abs(CL_-0.10)<0.005 );
              
               
    end
    theta_tune_is_train=(C_theta_1.* theta_tune_is_train);
end