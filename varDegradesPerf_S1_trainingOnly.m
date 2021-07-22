% this only works if everything is kept in memory from
% varDegradesPerf_S1_HallemOslenlnp.m

for randomTrials=1:ll
    randomTrials
    tic
    
    
    load(strcat('VarDegradesPerf_fly_wNoise_HOInp',num2str(randomTrials),num2str(noiseScale)));
    
    Activations=zeros(n,odors*numTrials);
    ActivationsHomog=zeros(n,odors*numTrials);
    Activations_varW_fixedN= zeros(n,odors*numTrials);
    Activations_varN_fixedW= zeros(n,odors*numTrials);
    Activations_Ftheta=zeros(n,odors*numTrials);
    ActivationsHomog_Ftheta=zeros(n,odors*numTrials);
    Activations_varW_fixedN_and_theta= zeros(n,odors*numTrials);
    Activations_varN_fixedW_and_theta= zeros(n,odors*numTrials);
    
    
    Y=zeros(n,odors*numTrials);
    YHomog=zeros(n,odors*numTrials);
    Y_varw_fixedN= zeros(n,odors*numTrials);
    Y_varN_fixedW= zeros(n,odors*numTrials);
    Y_Ftheta=zeros(n,odors*numTrials);
    YHomog_Ftheta=zeros(n,odors*numTrials);
    Y_varw_fixedN_and_theta= zeros(n,odors*numTrials);
    Y_varN_fixedW_and_theta= zeros(n,odors*numTrials);
    
    
    for trial = 1:(odors*numTrials)
        
        ActivationsHomog_Ftheta(:,trial) = thisW_HomogModel'*PNtrials(:,trial  );
        YHomog_Ftheta(:,trial)=(( ActivationsHomog_Ftheta(:,trial)-(APLgains(4))*repmat(sum(ActivationsHomog_Ftheta(:,trial),1),n,1)-thetaH_Ftheta)>0 ).*( ActivationsHomog_Ftheta(:,trial)-APLgains(4)*repmat(sum(ActivationsHomog_Ftheta(:,trial),1),n,1)-thetaH_Ftheta);
        
        Activations_Ftheta(:,trial) = thisW'*PNtrials(:,trial );
        Y_Ftheta(:,trial)=(( Activations_Ftheta(:,trial)-(APLgains(3))*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta)>0 ).*( Activations_Ftheta(:,trial)-APLgains(3)*repmat(sum(Activations_Ftheta(:,trial),1),n,1)-thetaS_Ftheta);
        
        Activations_varW_fixedN_and_theta(:,trial) = thisW_varWonly'*PNtrials(:,trial );
        Y_varw_fixedN_and_theta(:,trial)=(( Activations_varW_fixedN_and_theta(:,trial)-(APLgains(6) )*repmat(sum(Activations_varW_fixedN_and_theta(:,trial),1),n,1)-theta_varw_fixedN_scalar)>0 ).*( Activations_varW_fixedN_and_theta(:,trial)-APLgains(6)*repmat(sum(Activations_varW_fixedN_and_theta(:,trial),1),n,1)-theta_varw_fixedN_scalar);
        
        Activations_varN_fixedW_and_theta(:,trial) = thisW_varN'*PNtrials(:,trial );
        Y_varN_fixedW_and_theta(:,trial)=(( Activations_varN_fixedW_and_theta(:,trial)-(APLgains(8) )*repmat(sum(Activations_varN_fixedW_and_theta(:,trial),1),n,1)-theta_varN_fixedW_scalar)>0 ).*( Activations_varN_fixedW_and_theta(:,trial)-APLgains(8)*repmat(sum(Activations_varN_fixedW_and_theta(:,trial),1),n,1)-theta_varN_fixedW_scalar);
        
        ActivationsHomog(:,trial) = thisW_HomogModel'*PNtrials(:,trial );
        YHomog(:,trial)=(( ActivationsHomog(:,trial)-(APLgains(2))*repmat(sum(ActivationsHomog(:,trial),1),n,1)-thetaH)>0 ).*( ActivationsHomog(:,trial)-APLgains(2)*repmat(sum(ActivationsHomog(:,trial),1),n,1)-thetaH);
        
        Activations(:,trial) = thisW'*PNtrials(:,trial );
        Y(:,trial)=(( Activations(:,trial)-(APLgains(1))*repmat(sum(Activations(:,trial),1),n,1)-thetaS)>0 ).*( Activations(:,trial)-APLgains(1)*repmat(sum(Activations(:,trial),1),n,1)-thetaS);
        
        Activations_varW_fixedN(:,trial) = thisW_varWonly'*PNtrials(:,trial );
        Y_varw_fixedN(:,trial)=(( Activations_varW_fixedN(:,trial)-(APLgains(5) )*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN)>0 ).*( Activations_varW_fixedN(:,trial)-APLgains(5)*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN);
        
        Activations_varN_fixedW(:,trial) = thisW_varN'*PNtrials(:,trial );
        Y_varN_fixedW(:,trial)=(( Activations_varN_fixedW(:,trial)-(APLgains(7) )*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw)>0 ).*( Activations_varN_fixedW(:,trial)-APLgains(7)*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw);
        
    end
    
    logLRs = [-5 -4 -3 -2.75 -2.5 -2.25 -2 -1 0 1];
    for l_r=1:length(logLRs)
        
        % starting at initial random values for KC-MBon weights
        
        WopAllOdours=1*rand(n,2);
        WopAllOdoursHomog=WopAllOdours;
        WopAllOdours_Ftheta=WopAllOdours;
        WopAllOdoursHomog_Ftheta=WopAllOdours;
        
        WopAllOdoursVarW_fixedN= WopAllOdours;
        WopAllOdoursVarW_fixedN_and_theta=WopAllOdours;
        
        WopAllOdoursVarN_fixedW=WopAllOdours;
        WopAllOdoursVarN_fixedW_and_theta=WopAllOdours;
        
        WopFromPNs= 1*rand(m,2);
        %                 alpha= 0.001*(10^(l_r));
        alpha=10^(logLRs(l_r));
        
        c=1;
        ceq=1;
        ch=1;
        
        YHomogtemp=reshape(YHomog,n,odors,numTrials);
        YHomog_Fthetatemp=reshape(YHomog_Ftheta,n,odors,numTrials);
        Ytemp= reshape(Y,n,odors,numTrials);
        Y_Fthetatemp= reshape(Y_Ftheta,n,odors,numTrials);
        Y_varw_fixedNtemp= reshape(Y_varw_fixedN,n,odors,numTrials);
        Y_varw_fixedN_and_thetatemp= reshape(Y_varw_fixedN_and_theta,n,odors,numTrials);
        Y_varN_fixedWtemp= reshape(Y_varN_fixedW,n,odors,numTrials);
        Y_varN_fixedW_and_thetatemp=reshape(Y_varN_fixedW_and_theta,n,odors,numTrials);
        
        % rescale the KC responses to be from 0-1
        Ytemp=rescale(Ytemp);
        YHomogtemp=rescale(YHomogtemp);
        Y_Fthetatemp=rescale(Y_Fthetatemp);
        YHomog_Fthetatemp=rescale(YHomog_Fthetatemp);
        Y_varw_fixedNtemp=rescale(Y_varw_fixedNtemp);
        Y_varw_fixedN_and_thetatemp=rescale(Y_varw_fixedN_and_thetatemp);
        Y_varN_fixedWtemp=rescale(Y_varN_fixedWtemp);
        Y_varN_fixedW_and_thetatemp= rescale(Y_varN_fixedW_and_thetatemp);
        
        
        YHomogtr=YHomogtemp(:,:,1:numtrainingSamples);
        YHomog_Fthetatr=YHomog_Fthetatemp(:,:,1:numtrainingSamples);
        Ytr=Ytemp(:,:,1:numtrainingSamples);
        Y_Fthetatr=Y_Fthetatemp(:,:,1:numtrainingSamples);
        Y_varWtr= Y_varw_fixedNtemp(:,:,1:numtrainingSamples);
        Y_varWFtr=Y_varw_fixedN_and_thetatemp(:,:,1:numtrainingSamples);
        Y_varNtr= Y_varN_fixedWtemp(:,:,1:numtrainingSamples);
        Y_varNFtr= Y_varN_fixedW_and_thetatemp(:,:,1:numtrainingSamples);
        
        % learning the MBON- preceptron weights in the output layer in each of the 8 models: by synaptic depression
        % punishing the connection to the MBON decision (MBON avoidance & MBON approach) opposite to the odor class ('good'
        % rewarded odor and 'bad' punished odor)
        %  i.e.learning from mistakes
        
        
        for odour=1:odors
            
            if( ~ isempty(find(classAction1==odour)) )
                delta =  exp( -(alpha/mean(YHomogtr(:)))* sum(YHomogtr(:,odour,:),3) );
                deltaWH(:,ch)= WopAllOdoursHomog(:,2).*(delta-1);
                ch=ch+1;
                WopAllOdoursHomog(:,2)= WopAllOdoursHomog(:,2) .*delta;
                
            else
                delta = exp(- (alpha/mean(YHomogtr(:)))* sum(YHomogtr(:,odour,:),3) );
                WopAllOdoursHomog(:,1)= WopAllOdoursHomog(:,1) .*delta;
                
            end
            
            
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
            
            
            if( ~ isempty(find(classAction1==odour)) )
                
                delta =  exp( -(alpha/mean(Y_varWtr(:)))* sum(Y_varWtr(:,odour,:),3) );
                WopAllOdoursVarW_fixedN(:,2)= WopAllOdoursVarW_fixedN(:,2) .*delta;
                
            else
                delta = exp(- (alpha/mean(Y_varWtr(:)))* sum(Y_varWtr(:,odour,:),3) );
                WopAllOdoursVarW_fixedN(:,1)= WopAllOdoursVarW_fixedN(:,1) .*delta;
            end
            
            if( ~ isempty(find(classAction1==odour)) )
                
                delta =  exp( -(alpha/mean(Y_varWFtr(:)))* sum(Y_varWFtr(:,odour,:),3) );
                WopAllOdoursVarW_fixedN_and_theta(:,2)= WopAllOdoursVarW_fixedN_and_theta(:,2) .*delta;
                
            else
                delta = exp(- (alpha/mean(Y_varWFtr(:)))* sum(Y_varWFtr(:,odour,:),3) );
                WopAllOdoursVarW_fixedN_and_theta(:,1)= WopAllOdoursVarW_fixedN_and_theta(:,1) .*delta;
            end
            
            
            if( ~ isempty(find(classAction1==odour)) )
                
                delta =  exp( -(alpha/mean(Y_varNtr(:)))* sum(Y_varNtr(:,odour,:),3) );
                
                WopAllOdoursVarN_fixedW(:,2)= WopAllOdoursVarN_fixedW(:,2) .*delta;
                
            else
                delta = exp(- (alpha/mean(Y_varNtr(:)))* sum(Y_varNtr(:,odour,:),3) );
                WopAllOdoursVarN_fixedW(:,1)= WopAllOdoursVarN_fixedW(:,1) .*delta;
            end
            
            if( ~ isempty(find(classAction1==odour)) )
                
                delta =  exp( -(alpha/mean(Y_varNFtr(:)))* sum(Y_varNFtr(:,odour,:),3) );
                WopAllOdoursVarN_fixedW_and_theta(:,2)= WopAllOdoursVarN_fixedW_and_theta(:,2) .*delta;
                
            else
                delta = exp(- (alpha/mean(Y_varNFtr(:)))* sum(Y_varNFtr(:,odour,:),3) );
                WopAllOdoursVarN_fixedW_and_theta(:,1)= WopAllOdoursVarN_fixedW_and_theta(:,1) .*delta;
            end
            
            
        end
        
        %% perfromance as a function of the strictness of the decision making
        %% this strictness is dictated by C in the soft-max function.
        %% so given the same fly, same task, and after learning measure the performance as f(c)
        
        for c=1:Crange
            
            C=C_SoftMax*(10^c);
            [acc,accH,acc7, acc9]=KernelTesting(C,WopAllOdours,WopAllOdoursHomog, WopAllOdoursVarW_fixedN,WopAllOdoursVarN_fixedW , classAction1,numTrials,numtrainingSamples,Ytemp,YHomogtemp, Y_varw_fixedNtemp, Y_varN_fixedWtemp);
            
            test_p_ra(mods,randomTrials,noiseScale, l_r,c)=acc;
            test_p_raH(mods,randomTrials,noiseScale, l_r,c)=accH;
            test_p_ra_varW_FixedN(mods,randomTrials,noiseScale, l_r,c)= acc7;
            test_p_ra_varN_FixedW(mods,randomTrials,noiseScale, l_r,c)= acc9;
            
            [acc1,accH1,acc8,acc10]=KernelTesting_Ftheta(C,WopAllOdours_Ftheta,WopAllOdoursHomog_Ftheta,WopAllOdoursVarW_fixedN_and_theta, WopAllOdoursVarN_fixedW_and_theta,classAction1,numTrials,numtrainingSamples,Y_Fthetatemp,YHomog_Fthetatemp,Y_varw_fixedN_and_thetatemp,Y_varN_fixedW_and_thetatemp);
            
            test_p_ra_Fixedtheta(mods,randomTrials,noiseScale, l_r,c)=acc1;
            test_p_raH_FixedTheta(mods,randomTrials,noiseScale, l_r,c)=accH1;
            test_p_ra_varW_FixedN_FixedTheta(mods,randomTrials,noiseScale, l_r,c)= acc8;
            test_p_ra_varN_FixedW_FixedTheta(mods,randomTrials,noiseScale, l_r,c)=acc10;
            
            
        end
        
        
        
    end
    
    toc
    
end

save('test_p_ra.mat','test_p_ra');
save('test_p_raH.mat','test_p_raH');
save('test_p_ra_Fixedtheta.mat','test_p_ra_Fixedtheta');
save('test_p_raH_FixedTheta.mat','test_p_raH_FixedTheta');
save('test_p_ra_varW_FixedN_FixedTheta.mat','test_p_ra_varW_FixedN_FixedTheta')
save('test_p_ra_varW_FixedN.mat','test_p_ra_varW_FixedN')
save('test_p_ra_varN_FixedW_FixedTheta.mat','test_p_ra_varN_FixedW_FixedTheta')
save('test_p_ra_varN_FixedW.mat','test_p_ra_varN_FixedW')