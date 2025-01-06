function MBmodel = optimiseMBparams_homeostaticExcitation(MBmodel, PNactivity)
%% this is the "blue model"

DEBUG = true;

sigmoid_factor = 1.;
Sigmoid_derivative = @(x) exp(-x./sigmoid_factor)./((1.+exp(-x./sigmoid_factor)).^2) ;
Sigmoid = @(x) 1./(1.+exp(-x./sigmoid_factor)) ;

n = MBmodel.nKCs;
m = MBmodel.nPNs;
theta = MBmodel.theta;
weights = MBmodel.PNtoKC'; %take the transpose here to simplify
nResponses = size(PNactivity,2)*size(PNactivity,3);
PNtoKCmask = double(weights>0);

eta_Ctheta = 1; %scales adjustment steps for C_theta
C_theta=1; %scaling factor to achieve the correct average coding level (10 resp. 20%)

% eta_APLgain = 0.0000001; %scales adjustment steps for ALPgain
eta_APLgain = 0.0000001; %scales adjustment steps for ALPgain
APLgain = 0.000001;

lifetimeSparseness = 0.51;
sparseness_margin = lifetimeSparseness * 0.06;
eta_weights = 0.12;

constraints=0; %%
nLoops = 0;
% while the Spar.constraints are NOT TRUE: repeat
detectDeadEnd = true;
if detectDeadEnd
    APLgain_prev = APLgain;
    C_prev = C_theta;
    weights_prev = weights;
end

tic
while(~constraints)
    nLoops = nLoops+1;
    if nLoops>1e4
        error('Optimisation did not converge')
    end

    A = weights * PNactivity; %took transpose of PN-KC weights earlier
    Y_disInh = A - C_theta .* theta;
    CL_disInh = mean(Y_disInh(:)>0); %inhibition absent, should get to 20%
    
    %% we want to change the theta so to achieve CL_incInh=2xCL_disInh
	% c.f. eq.14 and previous
	% CL stands for Coding Level, and we want 10% resp. 20% of KCs active with/without inhibition
    dsig_dy = Sigmoid_derivative(Y_disInh);
    if any(isnan(dsig_dy(:)))
        fprintf('found NaNs in iteration %d',nLoops)
    end
    dsig_dy(isnan(dsig_dy)) = 0;
    depsi_dtheta= -(Y_disInh>0).* dsig_dy.* (repmat(theta,1,nResponses));
    Grad = (CL_disInh-0.20) .* mean(depsi_dtheta(:));
    C_theta = C_theta - eta_Ctheta.*Grad;
    
    if (C_theta<0)
        error('the scale factor in the random model is -ve')
    end
    
    %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
    % similar calculation, but including estimate of APLgain, and with adjusted C_theta
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    CL_incInh = mean(Y_incInh(:)>0);

	%c.f. eq.18 and previous, adjusting APLgain alpha
    dsig_dy = Sigmoid_derivative(Y_incInh);
    dsig_dy(isnan(dsig_dy)) = 0;
    dAct_dalpha = sum(A,1);
    dsig_dalpha = -(Sigmoid(Y_incInh)) .* (repmat(dAct_dalpha,n,1)) .*  dsig_dy;
	Grad_alpha = (CL_incInh-0.10) * mean(dsig_dalpha(:));
    APLgain = APLgain - eta_APLgain * Grad_alpha;

    %new part
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    Y_incInh(Y_incInh<0) = 0;
    avgAKcs = mean(Y_incInh,2); %lifetime average activity (mean across trials)
    errorInActivity = (avgAKcs-lifetimeSparseness); %what is lifetimeSparseness set to?
    weights = weights - eta_weights * PNtoKCmask .* repmat(errorInActivity,1,m);
    weights(weights<0) = 0;

    % sanity check
    if detectDeadEnd
        if (APLgain==APLgain_prev) & (C_theta==C_prev) & all(weights(:)==weights_prev(:))
            fprintf('loop %d',nLoops)
            error('Values remained unchanged without fulfilling the criteria!')
        end
        APLgain_prev = APLgain;
        C_prev = C_theta;
        weights_prev = weights;
    end

    %% checking if the constraints are met:
	% repeat calculations done during optimisation but with updated values C_theta, APLgain
    % first for zero-inhibition condition
    % skip calculating activation, it doesn't change
    A = weights * PNactivity; %took transpose of PN-KC weights earlier
    Y_disInh = A - C_theta .* theta ;
    CL_disInh = mean(Y_disInh(:)>0);
    
	%  then for KC activity with APL inhibition
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    CL_incInh = mean(Y_incInh(:)>0);
    Y_incInh(Y_incInh<0) = 0;
    avgAKcs = mean(Y_incInh,2); %additional condition
    %check conditions on CL with/without inhibition and lifetime sparseness
    constraints = all( abs(avgAKcs-lifetimeSparseness)<sparseness_margin ) & ( abs( (CL_disInh/CL_incInh) - 2.0)<0.2 ) &...
        (abs( CL_incInh-0.10)<0.01);

    if DEBUG
        disp([nLoops, CL_disInh, CL_incInh, mean(avgAKcs),std(avgAKcs), C_theta, APLgain])
    end

end

toc
disp(sprintf('Optimisation took %d loops',nLoops))

MBmodel.C_theta = C_theta;
MBmodel.alpha = APLgain;
MBmodel.PNtoKC = weights'; %undo the transpose iff it was done before loop
MBmodel.CL_disInh = CL_disInh;
MBmodel.CL_incInh = CL_incInh;
MBmodel.lifetimeActivity = avgAKcs;
MBmodel.modelType = 'homeostatic PN->KC weigths aka BLUE model';

end