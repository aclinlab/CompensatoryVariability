tic
A=zeros(n,odorsTuning_training*numtrainingSamples);
Y_d=zeros(n,odorsTuning_training*numtrainingSamples);

iterr=1;
codingLevelDummy=[];
conditions=0;
cw_=1;
C_ = 1.3;
APLgains_tune_is_train(4) = 6e-5;

if(tune~=3)
    A0=(0.4).*ones(n,1);
    epsilon= A0(1)*0.07;
    
else
    
    A0=(0.3).*ones(n,1);
    epsilon= A0(1)*0.07;
    
end
eta_0=0.05;
drop_1=0.7;
iterrDrop_1=1000;

eta_gradAct_theta_0=0.15;%0.05;
drop=0.7;
iterDrop=1000;
avgact_trace=[];

while(~conditions)
    
    % with inhibition gain absent, make sure that
    % CL=20% with the random values for theta.
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        A(:,trial) = (thisW_Kennedy)'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( A(:,trial)-(C_1.*theta_Activity_homeo_tune_is_train))>0 ).*( A(:,trial)-(C_1.*theta_Activity_homeo_tune_is_train));
        codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    
    InhAbs_CL=mean(codingLevelDummy);
    
    depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
    depsi1_dy(isnan(depsi1_dy))=0;
    
    depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_Activity_homeo_tune_is_train,1,odorsTuning_training*numtrainingSamples));
    
    %               *mean(Y_d(:));
    eta=eta_0*(drop_1^(floor(iterr/iterrDrop_1)));
    
    
    Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
    
    C_1=C_1-(eta_0*Grad);
    
    if (C_1<0)
        
        error('the scale factor in our magenta model, in tune is train is -ve!!')
    end
    
    
    %% CL=10% constraint
    eta_2=0.000000005;
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        Activations(:,trial) = (thisW_Kennedy)'*PNtrials_tune_train(:,trial);
        Y_(:,trial)=(( Activations(:,trial)-(APLgains_tune_is_train(4))*repmat(sum(Activations(:,trial),1),n,1)-(C_1.*theta_Activity_homeo_tune_is_train))>0 ).*( Activations(:,trial)-APLgains_tune_is_train(4)*repmat(sum(Activations(:,trial),1),n,1)-(C_1.*theta_Activity_homeo_tune_is_train));
        codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
    end
    CL_=mean(codingLevelDummy);
    
    dsig_dy=(exp(0.9.*Y_)./((1+exp(0.9.*Y_)).^2));
    dsig_dy(isnan(dsig_dy))=0;
    
    dAct_dalpha= sum(Activations,1);
    
    if ( any(isinf(dAct_dalpha)) )
        
        stopp=1;
    end
    
    dsig_dalpha= -(Y_>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
    
    
    Grad_alpha= ((CL_)-0.10)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(dsig_dalpha(:)));
    APLgains_tune_is_train(4)= APLgains_tune_is_train(4)- eta_2*(Grad_alpha);
    
    %% now tune theta for equal average activity level.
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        Activations(:,trial) = (thisW_Kennedy)'*PNtrials_tune_train(:,trial );
        Y_(:,trial)=(( Activations(:,trial)-(APLgains_tune_is_train(4))*repmat(sum(Activations(:,trial),1),n,1)-(C_1.*theta_Activity_homeo_tune_is_train))>0 ).*( Activations(:,trial)- APLgains_tune_is_train(4)*repmat(sum(Activations(:,trial),1),n,1)-(C_1.*theta_Activity_homeo_tune_is_train));
        
    end
    
    avgAKcs=mean(Y_,2);
    
    errorInActivity=(avgAKcs-A0);
    adEta=eta_gradAct_theta_0*(drop^(floor(iterr/iterDrop)));
    theta_Activity_homeo_tune_is_train= theta_Activity_homeo_tune_is_train-( ((eta_gradAct_theta_0).*((-1.*C_1.*(errorInActivity)))));
    theta_Activity_homeo_tune_is_train(theta_Activity_homeo_tune_is_train<0)=0;
    
    
    %% check the constraints
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        ActivationsDummy(:,trial) = (thisW_Kennedy)'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_1.*theta_Activity_homeo_tune_is_train))>0 ).*( ActivationsDummy(:,trial)-(C_1.*theta_Activity_homeo_tune_is_train));
        codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    
    InhAbs_CL=mean(codingLevelDummy);
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        Activations(:,trial) = (thisW_Kennedy)'*PNtrials_tune_train(:,trial);
        Y_(:,trial)=(( Activations(:,trial)-(APLgains_tune_is_train(4))*repmat(sum(Activations(:,trial),1),n,1)-(C_1.*theta_Activity_homeo_tune_is_train))>0 ).*( Activations(:,trial)-APLgains_tune_is_train(4)*repmat(sum(Activations(:,trial),1),n,1)-(C_1.*theta_Activity_homeo_tune_is_train));
        codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
        
    end
    CL_=mean(codingLevelDummy);
    avgAKcs=mean(Y_,2);
    
    
    TunedKCs=size(find(abs(avgAKcs-A0)<epsilon));
    
    conditions=  TunedKCs(1)>=1995 &( abs( (InhAbs_CL/CL_) - 2.0)<=0.3 ) &( (abs(CL_-0.10)) <=0.015 );
    if mod(iterr,10)==0
        disp('tune thetahomeo on training');
        disp(iterr)
        disp(CL_);
        disp(InhAbs_CL)
        disp(nnz(abs(avgAKcs-A0)<epsilon))
    end
    avgact_trace(:,iterr) = avgAKcs;

    iterr=iterr+1;
    
    
    
end

theta_Activity_homeo_tune_is_train=C_1.*theta_Activity_homeo_tune_is_train;
toc
disp('finished tune thetahomeo on training');