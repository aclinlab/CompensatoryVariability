function MBmodel = optimiseMBparams_APLgain_Ctheta(MBmodel, PNactivity)

% sigmoid_factor = -0.9;
% Sigmoid_derivative = @(x) exp(-sigmoid_factor.*x)./((1.+exp(-sigmoid_factor.*x)).^2) ;
% Sigmoid = @(x) 1./(1.+exp(-sigmoid_factor.*x)) ;
sigmoid_factor = -1.11;
Sigmoid_derivative = @(x) exp(-x./sigmoid_factor)./((1.+exp(-x./sigmoid_factor)).^2) ;
Sigmoid = @(x) 1./(1.+exp(-x./sigmoid_factor)) ;


%n=2000 KCs, numtrainingSamples=15
n = MBmodel.nKCs;
theta = MBmodel.theta;
thisW = MBmodel.PNtoKC;
nResponses = size(PNactivity,2)*size(PNactivity,3);

eta_Ctheta = 1; %scales adjustment steps for C_theta
C_theta=1; %scaling factor to achieve the correct average coding level (10 resp. 20%)

% eta_APLgain = 0.0000001; %scales adjustment steps for ALPgain
eta_APLgain = 0.0000001; %scales adjustment steps for ALPgain
APLgain = 0.000001;

conditions = false;
nLoops = 0;
detectDeadEnd = true;
if detectDeadEnd
    APLgain_prev = APLgain;
    C_prev = C_theta;
end
% while the Spar.constraints are NOT TRUE: repeat
while(~conditions)
    nLoops = nLoops+1;
    if nLoops>1e3
        error('Optimisation did not converge')
    end
    %calculate the KC responses to PN input for training trials
    A = thisW' * PNactivity;
    Y_disInh = A - C_theta .* theta;
    CL_disInh = mean(Y_disInh(:)>0); %inhibition absent, should get to 20%
    
    %% we want to change the theta so to achieve InhAbs_CL=2xCL
	% c.f. eq.14 and previous
	% CL stands for Coding Level, and we want 10% resp. 20% of KCs active with/without inhibition
    dsig_dy = Sigmoid_derivative(Y_disInh);
    if any(isnan(dsig_dy(:)))
        fprintf('found NaNs in iteration %d',nLoops)
    end
    dsig_dy(isnan(dsig_dy))=0;
    depsi_dtheta= -(Y_disInh>0).* dsig_dy.* (repmat(theta,1,nResponses));
    % Grad= ((CL_disInh)-0.20) ./ (n*nResponses) .* (sum(depsi_dtheta(:)));
    Grad = ((CL_disInh)-0.20) .* (mean(depsi_dtheta(:)));
    C_theta = C_theta - eta_Ctheta.*(Grad);
    if (C_theta<0)
        error('the scale factor in the random model is -ve')
    end
     
    %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
    % similar calculation, but including estimate of APLgain, and with adjusted C_theta
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    CL_incInh = mean(Y_incInh(:)>0);
    
	%c.f. eq.18 and previous, adjusting APLgain alpha
    dsig_dy = Sigmoid_derivative(Y_incInh);
    dsig_dy(isnan(dsig_dy))=0;
    dAct_dalpha= sum(A,1);
    % dsig_dalpha = -(Y_incInh>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
    dsig_dalpha = -(Sigmoid(Y_incInh)).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
    Grad_alpha = ((CL_incInh)-0.10)/(n*nResponses)*(sum(dsig_dalpha(:)));
    APLgain = APLgain - eta_APLgain*(Grad_alpha);
    
    % sanity check
    if detectDeadEnd
        if (APLgain==APLgain_prev) & (C_theta==C_prev)
            fprintf('loop %d',nLoops)
            error('Values remained unchanged without fulfilling the criteria!')
        end
        APLgain_prev = APLgain;
        C_prev = C_theta;
    end
	
    %% checking if the constraints are met:
	% repeat calculations done during optimisation but with updated values C_theta, APLgain
    %  first for zero-inhibition condition
    Y_disInh = A - C_theta .* theta;
    CL_disInh = mean(Y_disInh(:)>0);
	%  then for KC activity with APL inhibition
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    CL_incInh = mean(Y_incInh(:)>0); 
    conditions= ( abs( (CL_disInh/CL_incInh) - 2.0)<0.2 ) &( abs(CL_incInh-0.10)<0.01 );
    
end
disp(sprintf('Optimisation took %d loops',nLoops))
% end optimisation, store results
MBmodel.C_theta = C_theta;
MBmodel.alpha = APLgain;
MBmodel.CL_disInh = CL_disInh;
MBmodel.CL_incInh = CL_incInh;

end