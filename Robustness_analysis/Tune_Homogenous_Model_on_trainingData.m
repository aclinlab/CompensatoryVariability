function [thetaH_Ftheta_tune_train, APLgainP_tune_train_H,CL_, InhAbs_CL] = Tune_Homogenous_Model_on_trainingData(thisW_HomogModel,...
   PNtrials_tune_train, thetaH_Ftheta_tune_train, APLgainP_tune_train_H, odorsTuning_training,numtrainingSamples,C_thetaH_1)


n=2000;
m=24;
constraints=0;
A=zeros(n,odorsTuning_training*numtrainingSamples);
Y=zeros(n,odorsTuning_training*numtrainingSamples);
codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);

   while(~constraints)

                       eta=1;   
                       for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW_HomogModel'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)- (C_thetaH_1.*thetaH_Ftheta_tune_train))>0 ).*( A(:,trial)-(C_thetaH_1.*thetaH_Ftheta_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);

                       end
                       InhAbs_mHomog=mean(codingLevel);    
                       %% we want to change the theta so to achieve InhAbs_CL=2xCL


                       de1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                       de1_dy(isnan(de1_dy))=0;

                       de1_dtheta= -(Y>0).* de1_dy.*(repmat(thetaH_Ftheta_tune_train,n,odorsTuning_training*numtrainingSamples));


                       Grad= ((InhAbs_mHomog)-0.20)*(1/(n*odorsTuning_training*numtrainingSamples)).*(sum(de1_dtheta(:) ));

                       C_thetaH_1=C_thetaH_1 - eta.*(Grad);
                       
                       if (C_thetaH_1<0)
                           error('the scale factor in black model with tune is train, is -ve')
                       end


                        % allow alpha to be a free parameter
                       A=zeros(n,odorsTuning_training*numtrainingSamples);
                       Y=zeros(n,odorsTuning_training*numtrainingSamples);
                       codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);

                        eta_2=0.0000001;

                       for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW_HomogModel'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(APLgainP_tune_train_H)*repmat(sum(A(:,trial),1),n,1)-(C_thetaH_1.*thetaH_Ftheta_tune_train))>0 ).*( A(:,trial)-APLgainP_tune_train_H*repmat(sum(A(:,trial),1),n,1)-(C_thetaH_1.*thetaH_Ftheta_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);

                       end
                       CLAllFixed=mean(codingLevel);

                       dsig_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                       dsig_dy(isnan(dsig_dy))=0;

                       dAct_dalpha= sum(A,1);

                       dsig_dalpha= -(Y>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;


                       Grad_alpha= ((CLAllFixed)-0.10)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(dsig_dalpha(:)));
                       APLgainP_tune_train_H= APLgainP_tune_train_H- eta_2*(Grad_alpha);



                        %% check if the constraints are satisfied 
                       A=zeros(n,odorsTuning_training*numtrainingSamples);
                       Y=zeros(n,odorsTuning_training*numtrainingSamples);
                       codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);


                        for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW_HomogModel'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_thetaH_1.*thetaH_Ftheta_tune_train))>0 ).*( A(:,trial)-(C_thetaH_1.*thetaH_Ftheta_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n); 
                        end

                         InhAbs_CL=mean(codingLevel);
                       %
                       
                       A=zeros(n,odorsTuning_training*numtrainingSamples);
                       Y=zeros(n,odorsTuning_training*numtrainingSamples);
                       codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);

                          for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW_HomogModel'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(APLgainP_tune_train_H)*repmat(sum(A(:,trial),1),n,1)-(C_thetaH_1.*thetaH_Ftheta_tune_train))>0 ).*( A(:,trial)-APLgainP_tune_train_H*repmat(sum(A(:,trial),1),n,1)-(C_thetaH_1.*thetaH_Ftheta_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n); 

                          end
                         CL_=mean(codingLevel);


                         constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.1) &( abs(CL_-0.10)<0.01 );
              

   end
    thetaH_Ftheta_tune_train=(C_thetaH_1.*thetaH_Ftheta_tune_train);

end