Activations =zeros(n,odors*numtrainingSamples);
ActivationsDummy= zeros(n, odors*numtrainingSamples);
Conn=zeros(m,n);
Conn(find(thisW_ActivityBasedComp_noxjk_tuneis_train))=1;
mask= zeros (m,n);
mask(find(thisW_ActivityBasedComp_noxjk_tuneis_train))=1;
C_=1;
APLgains_noxjk_tuneis_train=0;
Inhabs_CLV=[];
CLV=[];
conditions=0; %%
A=zeros(n,odors*numtrainingSamples);
Y_d=zeros(n,odors*numtrainingSamples);
Y_=[];
codingLevelDummy=[];
tic
t=1;
while(~conditions)
    
    
    for trial = 1:(odors*numtrainingSamples)
        
        A(:,trial) = thisW_ActivityBasedComp_noxjk_tuneis_train'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( A(:,trial)-(C_.*theta_comp2_noxjk_tuneis_train) )>0 ).*( A(:,trial)-(C_.*theta_comp2_noxjk_tuneis_train));
        codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    InhAbs_CL=mean(codingLevelDummy);
    depsi1_dy=(exp(0.9.*Y_d)./((1+exp(0.9.*Y_d)).^2));
    depsi1_dy(isnan(depsi1_dy))=0;
    depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_comp2_noxjk_tuneis_train,1,odors*numtrainingSamples));
    eta=10;
    Grad= ((InhAbs_CL)-0.20)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
    C_=C_ - (eta*Grad);
    
    if (C_<0)
        error('the scale factor in the blue model is -ve!!')
    end
    
    
    eta_2=0.00000005;
    
    for trial = 1:(odors*numtrainingSamples)
        
        Activations(:,trial) = thisW_ActivityBasedComp_noxjk_tuneis_train'*PNtrials_tune_train(:,trial);
        Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk_tuneis_train)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk_tuneis_train))>0 ).*( Activations(:,trial)-APLgains_noxjk_tuneis_train*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk_tuneis_train));
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
    Grad_alpha= ((CL_)-0.10)*(1/(n*odors*numtrainingSamples))*(sum(dsig_dalpha(:)));
    APLgains_noxjk_tuneis_train= APLgains_noxjk_tuneis_train- eta_2*(Grad_alpha);
    
    
    %% tuning KC input weights to achieve target average activity A0
    
    for trial = 1:(odors*numtrainingSamples)
        Activations(:,trial) = thisW_ActivityBasedComp_noxjk_tuneis_train'*PNtrials_tune_train(:,trial );
        Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk_tuneis_train)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk_tuneis_train))>0 ).*( Activations(:,trial)-APLgains_noxjk_tuneis_train*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk_tuneis_train));
    end
    avgAKcs=mean(Y_,2);
    errorInActivity=(1).*repmat((avgAKcs-A0)',m,1);
    thisW_ActivityBasedComp_noxjk_tuneis_train= thisW_ActivityBasedComp_noxjk_tuneis_train-(0.05).*((1.*(mask.*errorInActivity)));
    
    %catch the -ve weights values
    if (~isempty(find(isinf(thisW_ActivityBasedComp_noxjk_tuneis_train) )))
        g=1;
    end
    thisW_ActivityBasedComp_noxjk_tuneis_train(find(thisW_ActivityBasedComp_noxjk_tuneis_train<0))=0;
    
    % check if constraints are satisfied
    for trial = 1:(odors*numtrainingSamples)
        
        ActivationsDummy(:,trial) = thisW_ActivityBasedComp_noxjk_tuneis_train'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( ActivationsDummy(:,trial)-(C_.*theta_comp2_noxjk_tuneis_train))>0 ).*( ActivationsDummy(:,trial)-(C_.*theta_comp2_noxjk_tuneis_train));
        codingLevelDummy(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    InhAbs_CL=mean(codingLevelDummy);
    
    for trial = 1:(odors*numtrainingSamples)
        Activations(:,trial) = thisW_ActivityBasedComp_noxjk_tuneis_train'*PNtrials_tune_train(:,trial);
        Y_(:,trial)=(( Activations(:,trial)-(APLgains_noxjk_tuneis_train)*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk_tuneis_train))>0 ).*( Activations(:,trial)-APLgains_noxjk_tuneis_train*repmat(sum(Activations(:,trial),1),n,1)-(C_.*theta_comp2_noxjk_tuneis_train));
        codingLevelDummy(trial)=  (sum(Y_(:,trial)>0,1)/n);
    end
    CL_=mean(codingLevelDummy);
    
    TunedKCs= size(find((abs(avgAKcs-A0)<epsilon)));
    conditions= TunedKCs(1)>=1995 &( abs( (InhAbs_CL/CL_) - 2.0)<0.2 ) &( (abs(CL_-0.10)) <=0.01 );
    
    if mod(t,10)==0
        disp('tune noxjk on training');
        disp(t)
        disp(CL_);
        disp(nnz(abs(avgAKcs-A0)<epsilon))
    end
    t=t+1;
    
end

theta_comp2_noxjk_tuneis_train=(C_.*theta_comp2_noxjk_tuneis_train);
toc
disp('finished tune noxjk on training');