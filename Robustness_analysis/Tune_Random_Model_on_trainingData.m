function [thetaS_tune_train, APLgainP_tune_train_S,CL_, InhAbs_CL] = Tune_Random_Model_on_trainingData(thisW,...
                             PNtrials_tune_train, thetaS_tune_train, APLgainP_tune_train_S,odorsTuning_training,numtrainingSamples,C_thetaS_1)
                         

n=2000;
m=24;                        
constraints=0; 

A=zeros(n,odorsTuning_training*numtrainingSamples);
Y=zeros(n,odorsTuning_training*numtrainingSamples);
eta=1; 
codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);

    while(~constraints)


                       for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_thetaS_1.*thetaS_tune_train))>0 ).*( A(:,trial)-(C_thetaS_1.*thetaS_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);

                       end
                       InhAbs_mSimp_=mean(codingLevel);    
                       
                       %% we want to change the theta so to achieve InhAbs_CL=2xCL

                       de1_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                       de1_dy(isnan(de1_dy))=0;

                       de1_dtheta= -(Y>0).* de1_dy.* (repmat(thetaS_tune_train,1,odorsTuning_training*numtrainingSamples));


                       Grad= ((InhAbs_mSimp_)-0.20)*(1/(n*odorsTuning_training*numtrainingSamples)).*(sum(de1_dtheta(:)));

                       C_thetaS_1=C_thetaS_1 - eta.*(Grad);

                      
                      if (C_thetaS_1<0)
                         error('the scale factor in the random model tune is train is -ve')
                     end

                 
                        %% allow alpha to be a free parameter
                        
                        eta_2=0.0000001;
                        A=zeros(n,odorsTuning_training*numtrainingSamples);
                        Y=zeros(n,odorsTuning_training*numtrainingSamples);
                        codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);

                       for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(APLgainP_tune_train_S)*repmat(sum(A(:,trial),1),n,1)-(C_thetaS_1.* thetaS_tune_train))>0 ).*( A(:,trial)-APLgainP_tune_train_S*repmat(sum(A(:,trial),1),n,1)-(C_thetaS_1.*thetaS_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n);

                       end
                       CLRand=mean(codingLevel);

                       ds_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2));
                       ds_dy(isnan(ds_dy))=0;

                       dAct_dalpha= sum(A,1);

                       ds_dalpha= -(Y>0).*(repmat(dAct_dalpha,n,1)).*  ds_dy;


                       Grad_alpha= ((CLRand)-0.10)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(ds_dalpha(:)));
                       APLgainP_tune_train_S= APLgainP_tune_train_S- eta_2*(Grad_alpha);


                       % check constraints
                        A=zeros(n,odorsTuning_training*numtrainingSamples);
                        Y=zeros(n,odorsTuning_training*numtrainingSamples);
                        codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);


                        for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(C_thetaS_1.* thetaS_tune_train))>0 ).*( A(:,trial)-(C_thetaS_1.* thetaS_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n); 
                        end

                         InhAbs_CL=mean(codingLevel);
                         
                        A=zeros(n,odorsTuning_training*numtrainingSamples);
                        Y=zeros(n,odorsTuning_training*numtrainingSamples);
                        codingLevel=zeros(1,odorsTuning_training*numtrainingSamples);

                          for trial = 1:(odorsTuning_training*numtrainingSamples) 

                            A(:,trial) = thisW'*PNtrials_tune_train(:,trial);
                            Y(:,trial)=(( A(:,trial)-(APLgainP_tune_train_S)*repmat(sum(A(:,trial),1),n,1)-(C_thetaS_1.* thetaS_tune_train))>0 ).*( A(:,trial)-APLgainP_tune_train_S*repmat(sum(A(:,trial),1),n,1)-(C_thetaS_1.* thetaS_tune_train));
                            codingLevel(trial)=  (sum(Y(:,trial)>0,1)/n); 

                          end
                         CL_=mean(codingLevel);


                         constraints= ( abs( (InhAbs_CL/CL_) - 2.0)<0.1 ) &( abs(CL_-0.10)<0.006 );

                          if(abs(InhAbs_CL-0.20)<0.000001)
                              break;
                          end


    end
    thetaS_tune_train=(C_thetaS_1.* thetaS_tune_train);
end