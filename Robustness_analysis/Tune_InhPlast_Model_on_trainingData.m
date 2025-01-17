tic
if tune==2||tune==1
    A0=(0.45).*ones(n,1);
    epsilon= A0(1)*0.07;
elseif(tune~=3)
    A0=(0.5).*ones(n,1);
    epsilon= A0(1)*0.07;
else
    A0=(0.4).*ones(n,1);
    epsilon= A0(1)*0.07;
end

conditions=0;

A=zeros(n,odorsTuning_training*numtrainingSamples);
Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
codingLevel=[];

iterr=1;
avgActTrace=[];
iter_till_exit=0;
while(~conditions && (iter_till_exit<10000))
    
    iter_till_exit=iter_till_exit+1;
    % with inhibition gain absent
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        A(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( A(:,trial)-(C_1.*theta_inhibitionPlast_tune_is_train))>0 ).*( A(:,trial)-(C_1.*theta_inhibitionPlast_tune_is_train));
        codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    
    InhAbs_CL=mean(codingLevel);
    
    depsi1_dy=(exp(0.91.*Y_d)./((1+exp(0.91.*Y_d)).^2));
    depsi1_dy(isnan(depsi1_dy))=0;
    
    depsi1_dtheta= -(Y_d>0).* depsi1_dy.* (repmat(theta_inhibitionPlast_tune_is_train,1,odorsTuning_training*numtrainingSamples));
    
    if tune==4
        eta=50;
    elseif tune==2 || tune==1
        eta=50;
    else
    eta=15; %*(0.7^(floor(iterr/2000)));
    end
    
    Grad= ((InhAbs_CL)-0.20)*(1/(n*odorsTuning_training*numtrainingSamples))*(sum(depsi1_dtheta(:) ));
    
    C_1=C_1 -(eta.* mean(Y_d(:)).*Grad);
    
    if (C_1<0)
        error('the scale factor in the green model, tune is train is -ve!')
        
    end
    
    %
    %A=zeros(n,odorsTuning_training*numtrainingSamples);
    Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
    codingLevel=[];
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        %A(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( A(:,trial)-(APLgains_model6_tune_is_train.*(sum(A(:,trial),1)))-(C_1.*theta_inhibitionPlast_tune_is_train))>0 ).*( A(:,trial)-(APLgains_model6_tune_is_train.*(sum(A(:,trial),1)))-(C_1.*theta_inhibitionPlast_tune_is_train));
        codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    CL_=mean(codingLevel);
    
    avgAKcs=mean(Y_d,2);
    errorInActivity=avgAKcs-A0;
    
    dsig_dy=(exp(0.91.*Y_d)./((1+exp(0.91.*Y_d)).^2));
    dsig_dy(isnan(dsig_dy))=0;
    
    xk_=repmat(sum(A),n,1);
    
    D_yj_alphaj= (-1.*mean(xk_,2));
    
    
    dAct_dalpha= sum(A,1);
    
    dsig_dalpha= -(Y_d>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
    
    dYik_dalphai= -(repmat(dAct_dalpha,n,1));
    
    if tune==4
        eta_o1 = 0.002;
    else
    eta_o1= 0.01; %0.15; %0.002;
    end
    eta_o2= 5e-7; %0.00000005;
    
    Grad_alpha1= ( ((eta_o1)) .*((CL_)-0.10)*(1/(n*odorsTuning_training*numtrainingSamples)).*(sum(dsig_dalpha,2)) );
    
    Grad_alpha2=( eta_o2.* (errorInActivity).*(D_yj_alphaj) );
    
    Grad_alpha= Grad_alpha1+Grad_alpha2 ;
    
    %                          eta_gradActiv_alpha=eta_gradAct_alpha_0*(drop^(floor(iterr/iterDrop)));
    
    APLgains_model6_tune_is_train= APLgains_model6_tune_is_train- 0.001.*((Grad_alpha+ ((0.000000001)*randn(1)) ));
    
    
    %% check constraints
    
    A=zeros(n,odorsTuning_training*numtrainingSamples);
    Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
    codingLevel=[];
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        A(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( A(:,trial)-(C_1.*theta_inhibitionPlast_tune_is_train))>0 ).*( A(:,trial)-(C_1.*theta_inhibitionPlast_tune_is_train));
        codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n);
    end
    
    InhAbs_CL=mean(codingLevel);
    
    A=zeros(n,odorsTuning_training*numtrainingSamples);
    Y_d=zeros(n,odorsTuning_training*numtrainingSamples);
    codingLevel=[];
    
    for trial = 1:(odorsTuning_training*numtrainingSamples)
        
        A(:,trial) = thisW_ActivityBasedComp_inhibitionPlast'*PNtrials_tune_train(:,trial);
        Y_d(:,trial)=(( A(:,trial)-(APLgains_model6_tune_is_train.*sum(A(:,trial),1))-(C_1.*theta_inhibitionPlast_tune_is_train))>0 ).*( A(:,trial)-(APLgains_model6_tune_is_train.*sum(A(:,trial),1))-(C_1.*theta_inhibitionPlast_tune_is_train));
        codingLevel(trial)=  (sum(Y_d(:,trial)>0,1)/n);
        
    end
    CL_=mean(codingLevel);
    avgAKcs=mean(Y_d,2);
    
    TunedKCs= size(find((abs(avgAKcs-A0)<epsilon)));
    %
    conditions=(TunedKCs(1)>=1995)  &( (abs(round(CL_,3)-0.10)) <=0.01) & ( round( abs( ((InhAbs_CL/CL_)) - 2.0),1) <=0.2 );
    
    if mod(iterr,10)==0
        disp('tune inhplast on training');
        disp(iterr)
        disp(CL_);
        disp(InhAbs_CL)
        disp(nnz(abs(avgAKcs-A0)<epsilon))
    end
    avgActTrace(:,iterr)=avgAKcs;
    
    iterr=iterr+1;
    
end
