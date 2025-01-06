function KCresponse = calculateKCresponse(MBmodel, PNactivity)
%% refers to the spiking/output activity of KCs, incorporating APL inhibition and spiking thresholds
% expects MBmodel to have the weighted connectivity matrix PNtoKC, alpha
% (APL gain) and C_theta (scaling for thresholds). PNactivity expected to
% be M-by-K, with M number of PNs for K odors

% excitInput variable refers to the first term in eq.3, i.e. the excitatory 
% drive to KCs before inhibition or thresholding
excitInput = MBmodel.PNtoKC' * PNactivity;
totalExc = sum(excitInput,1);
KCresponse = excitInput - MBmodel.alpha*totalExc - MBmodel.C_theta.*MBmodel.theta;
KCresponse(KCresponse<0) = 0;

% Activation variable refers to the first term in eq.3, i.e. the excitatory drive to KCs before inhibition or thresholding
% This is done first because the spiking activity combines the exc. and the APL inhib., whihc depends on total exc. as well
% Activations=zeros(MBmodel_varN_varW_varTh.nKCs,odors*numTrials);
% Y=zeros(MBmodel.nKCs,odors*numTrials); %model 1
% for trial = 1:size(PNactivity,2)
%     Activations(:,trial) = MBmodel.PNtoKC'*PNactivity(:,trial );
%     Y(:,trial)=(( Activations(:,trial)-(MBmodel.APLgain)*repmat(sum(Activations(:,trial),1),MBmodel.nKCs,1)-MBmodel.C_theta.*MBmodel.theta)>0 ).*( Activations(:,trial)-MBmodel.APLgain*repmat(sum(Activations(:,trial),1),MBmodel.nKCs,1)-MBmodel.theta);
%     % Y_varw_fixedN(:,trial)=(( Activations_varW_fixedN(:,trial)-(APLgains(5) )*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN)>0 ).*( Activations_varW_fixedN(:,trial)-APLgains(5)*repmat(sum(Activations_varW_fixedN(:,trial),1),n,1)-theta_varw_fixedN);
% 
%     % Activations_varN_fixedW(:,trial) = thisW_varN'*PNtrials(:,trial );
%     % Y_varN_fixedW(:,trial)=(( Activations_varN_fixedW(:,trial)-(APLgains(7) )*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw)>0 ).*( Activations_varN_fixedW(:,trial)-APLgains(7)*repmat(sum(Activations_varN_fixedW(:,trial),1),n,1)-theta_varN_fixedw);
% 
% end
end