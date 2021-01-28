tic

A=zeros(n,odorsTuning_training*numtrainingSamples);
Y_d=zeros(n,odorsTuning_training*numtrainingSamples);

codingLevel=[];

Conn=zeros(m,n);
Conn(find(thisW_ActivityBasedComp_tune_is_train))=1;

A0=(0.51)*ones(n,1);
epsilon= A0(1)*0.07; 
conditions=0;

                     
while(~conditions)
                   

                   eta=10;
                  % with inhibition gain absent
                   for trial = 1:(odorsTuning_training*numtrainingSamples) 
                   
                    A(:,trial) = thisW_ActivityBasedComp_tune_is_train'*PNtrials_tune_train(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_1.*theta_comp2_tune_is_train))>0 ).*( A(:,trial)-(C_1.*theta_comp2_tune_is_train));
                    codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n); 
                   end
                   
                   InhAbs_CL=mean(codingLevel);
                   
                  de1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                  de1_dy(isnan(de1_dy))=0;
                  
                  de1_dtheta= -(Y_d>0).* de1_dy.* (repmat(theta_comp2_tune_is_train,1,odorsTuning_training*numtrainingSamples));
                  
                  
                  Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(de1_dtheta(:) ));
                   
                  C_1=C_1 - eta*(Grad);
                  
                 if (C_1<0)
                     error('the scale factor in the blue model tune is train, is -ve!!')
                 end
                  
                   eta_2=0.00000005;
                  A=zeros(n,odorsTuning_training*numtrainingSamples);
                  Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
                  codingLevel=[];
                  
                  for trial = 1:(odorsTuning_training*numtrainingSamples) 
                   
                    A(:,trial) = thisW_ActivityBasedComp_tune_is_train'*PNtrials_tune_train(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(APLgains_tune_is_train(3))*repmat(sum(A(:,trial),1),n,1)-(C_1.*theta_comp2_tune_is_train))>0 ).*( A(:,trial)-APLgains_tune_is_train(3)*repmat(sum(A(:,trial),1),n,1)-(C_1.*theta_comp2_tune_is_train));
                    codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n); 

                  end
                 CL_=mean(codingLevel);
                 
                 dsig_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
                 dsig_dy(isnan(dsig_dy))=0;
                 dAct_dalpha= sum(A,1);
                   
                  if ( any(isinf(dAct_dalpha)) )
                      
                      stopp=1;
                  end
                  
                  dsig_dalpha= -(Y_d>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                  

                  Grad_alpha= ((CL_)-0.10)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(dsig_dalpha(:)));
                  APLgains_tune_is_train(3)= APLgains_tune_is_train(3)- eta_2*(Grad_alpha);
                 
                  
                  %% now, third optimization step, is to find (given this theta and alpha) 
                  %% the optimum weights to minimize the activity equalization error
                  A=zeros(n,odorsTuning_training*numtrainingSamples);
                  Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
                  codingLevel=[];
                  
                  
                  for trial = 1:(odorsTuning_training*numtrainingSamples)
                    A(:,trial) = thisW_ActivityBasedComp_tune_is_train'*PNtrials_tune_train(:,trial );
                    Y_d(:,trial)=(( A(:,trial)-(APLgains_tune_is_train(3))*repmat(sum(A(:,trial),1),n,1)-(C_1.*theta_comp2_tune_is_train))>0 ).*( A(:,trial)-APLgains_tune_is_train(3)*repmat(sum(A(:,trial),1),n,1)-(C_1.*theta_comp2_tune_is_train));

                  end


                    
                   avgAKcs=mean(Y_d,2);
                   errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);


                    Conn2=repmat(Conn,1,1,odorsTuning_training*numtrainingSamples);
                    Xik=reshape(PNtrials_tune_train(:,1:odorsTuning_training*numtrainingSamples),m,1,odorsTuning_training*numtrainingSamples);
                    Xik_=repmat(Xik,1,n,1);

                    filt_=  ( ((1-(InhibitionGain/n)) .* (Xik_.* Conn2)) ) ;

                    filt_Xjm= mean(filt_,3);
                    mask=filt_Xjm;

                   %% target activity is to change weight for already firing cells
                   %% if a KC is silent, leave it as is. and don't equalize its' excitability.
                   
                  
                  thisW_ActivityBasedComp_tune_is_train= thisW_ActivityBasedComp_tune_is_train-(0.05).*((1.*(mask.*errorInActivity)));                 
                 
                  %catch the -ve weights values
                  if (~isempty(find(isinf(thisW_ActivityBasedComp_tune_is_train) )))
                      g=1;
                  end
                  thisW_ActivityBasedComp_tune_is_train(find(thisW_ActivityBasedComp_tune_is_train<0))=0;
                 
                  
                  
                  A=zeros(n,odorsTuning_training*numtrainingSamples);
                  Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
                  codingLevel=[];
                  
                  for trial = 1:(odorsTuning_training*numtrainingSamples) 
                   
                    A(:,trial) = thisW_ActivityBasedComp_tune_is_train'*PNtrials_tune_train(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(C_1.*theta_comp2_tune_is_train))>0 ).*( A(:,trial)-(C_1.*theta_comp2_tune_is_train));
                    codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n); 
                  end
                 
                 InhAbs_CL=mean(codingLevel);

                 
                  A=zeros(n,odorsTuning_training*numtrainingSamples);
                  Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
                  codingLevel=[];
                  
                  for trial = 1:(odorsTuning_training*numtrainingSamples) 
                   
                    A(:,trial) = thisW_ActivityBasedComp_tune_is_train'*PNtrials_tune_train(:,trial);
                    Y_d(:,trial)=(( A(:,trial)-(APLgains_tune_is_train(3))*repmat(sum(A(:,trial),1),n,1)-(C_1.*theta_comp2_tune_is_train))>0 ).*( A(:,trial)-APLgains_tune_is_train(3)*repmat(sum(A(:,trial),1),n,1)-(C_1.*theta_comp2_tune_is_train));
                    codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n); 

                  end
                 CL_=mean(codingLevel);
                 avgAKcs= mean(Y_d,2);
                 
                
                 conditions= all(abs(avgAKcs-A0)<epsilon) &( abs( (InhAbs_CL/CL_) - 2.0)<0.1 ) &( (abs(CL_-0.10)) <=0.005 );
              
                

end

 theta_comp2_tune_is_train=(C_1.*theta_comp2_tune_is_train);                      