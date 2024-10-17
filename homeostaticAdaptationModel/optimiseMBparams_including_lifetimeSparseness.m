function MBmodel = optimiseMBparams_including_lifetimeSparseness(MBmodel, PNactivity)
%% this is the "magenta model"

% sigmoid_factor = -0.9;
% Sigmoid_derivative = @(x) exp(-sigmoid_factor.*x)./((1.+exp(-sigmoid_factor.*x)).^2) ;
% Sigmoid = @(x) 1./(1.+exp(-sigmoid_factor.*x)) ;
sigmoid_factor = -1.11;
Sigmoid_derivative = @(x) exp(-x./sigmoid_factor)./((1.+exp(-x./sigmoid_factor)).^2) ;
Sigmoid = @(x) 1./(1.+exp(-x./sigmoid_factor)) ;

%n=2000 KCs, numtrainingSamples=15
n = MBmodel.nKCs;
theta = MBmodel.theta;
weights = MBmodel.PNtoKC;
nResponses = size(PNactivity,2)*size(PNactivity,3);

% A = zeros(n,nResponses); %A is excitation to KC_i from PNs to odor k
% Y = zeros(n,nResponses); %resulting activity of KC without APL feedback
eta_Ctheta = 1; %scales adjustment steps for C_theta
C_theta=1; %scaling factor to achieve the correct average coding level (10 resp. 20%)

A0 = (0.51) .* ones(n,1);
epsilon = A0(1) * 0.06;

% eta_APLgain = 0.0000001; %scales adjustment steps for ALPgain
eta_APLgain = 0.0000001; %scales adjustment steps for ALPgain
APLgain = 0.000001;

constraints=0; %%
nLoops = 0;
% while the Spar.constraints are NOT TRUE: repeat
detectDeadEnd = true;
if detectDeadEnd
    APLgain_prev = APLgain;
    C_prev = C_theta;
    theta_prev = theta;
end

tic
while(~constraints)
    nLoops = nLoops+1;
    if nLoops>1e3
        error('Optimisation did not converge')
    end
    %calculate the KC responses to PN input for training trials
    A = weights' * PNactivity;
    Y = A - C_theta .* theta;
    CL_disInh = mean(Y(:)>0); %inhibition absent, should get to 20%
    
    %% we want to change the theta so to achieve CL_incInh=2xCL_disInh
	% c.f. eq.14 and previous
	% CL stands for Coding Level, and we want 10% resp. 20% of KCs active with/without inhibition
    % dsig_dy=(exp(0.9.*Y)./((1+exp(0.9.*Y)).^2)); %shouldn't this be 0.8 instead?
    dsig_dy = Sigmoid_derivative(Y);
    if any(isnan(dsig_dy(:)))
        fprintf('found NaNs in iteration %d',nLoops)
    end
    dsig_dy(isnan(dsig_dy))=0;
    depsi_dtheta= -(Y>0).* dsig_dy.* (repmat(theta,1,nResponses));
    % Grad= ((CL_disInh)-0.20) ./ (n*nResponses) .* (sum(depsi_dtheta(:)));
    Grad= ((CL_disInh)-0.20) .* (mean(depsi_dtheta(:)));
    C_theta=C_theta - eta_Ctheta.*(Grad);
    
    if (C_theta<0)
        error('the scale factor in the random model is -ve')
    end
    
    %% allow alpha to be a free parameter: achieving the 2nd spar. constraint CL=10%
    % similar calculation, but including estimate of APLgain, and with adjusted C_theta
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    CL_incInh = mean(Y_incInh(:)>0);
    
	%c.f. eq.18 and previous, adjusting APLgain alpha
    % dsig_dy=(exp(0.9.*YEqualizeddummy)./((1+exp(0.9.*YEqualizeddummy)).^2));
    dsig_dy = Sigmoid_derivative(Y_incInh);
    dsig_dy(isnan(dsig_dy))=0;
    dAct_dalpha= sum(A,1);
    % dsig_dalpha = -(Y_incInh>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
    dsig_dalpha = -(Sigmoid(Y_incInh)) .* (repmat(dAct_dalpha,n,1)) .*  dsig_dy;
    Grad_alpha = ((CL_incInh)-0.10)/(n*nResponses)*(sum(dsig_dalpha(:)));
    APLgain = APLgain - eta_APLgain*(Grad_alpha);
    
    %new part
    Y_incInh = A - APLgain.*repmat(sum(A,1),n,1) - C_theta.*repmat(theta,1,nResponses);
    avgAKcs = mean(Y_incInh,2); %lifetime average activity (mean actoss trials)
    errorInActivity = (avgAKcs-A0); %what is A0 set to?
    theta = theta - 0.01 * (-1.*C_theta .* errorInActivity);
    theta(theta<0) = 0;

    % sanity check
    if detectDeadEnd
        if (APLgain==APLgain_prev) & (C_theta==C_prev) & all(theta==theta_prev)
            fprintf('loop %d',nLoops)
            error('Values remained unchanged without fulfilling the criteria!')
        end
        APLgain_prev = APLgain;
        C_prev = C_theta;
        theta_prev = theta;
    end

    %% checking if the constraints are met:
	% repeat calculations done during optimisation but with updated values C_theta, APLgain
    % first for zero-inhibition condition
    % skip calculating activation, it doesn't change
    Y_disInh = A - C_theta .* theta ;
    CL_disInh = mean(Y_disInh(:)>0);
    
	%  then for KC activity with APL inhibition
    Y_incInh = A - APLgain .* sum(A,1) - C_theta .* theta;
    CL_incInh = mean(Y_incInh(:)>0);

    avgAKcs = mean(Y_incInh,2); %additional condition
    %check conditions on CL with/without inhibition and lifetime sparseness
    constraints = all( abs(avgAKcs-A0)<epsilon ) & ( abs( (CL_disInh/CL_incInh) - 2.0)<0.2 ) &...
        (abs( CL_incInh-0.10)<0.01);

    %disp([InhAbs_CL CL_ C_theta  APLgain])
    
end
toc
disp(sprintf('Optimisation took %d loops',nLoops))

MBmodel.C_theta = C_theta;
MBmodel.alpha = APLgain;
MBmodel.theta = theta;
MBmodel.CL_disInh = CL_disInh;
MBmodel.CL_incInh = CL_incInh;

% end optimisation, store results
% CLevelP(end+1)=CL_; %appends to vector, empty until now (renewed for each noise level)
% theta=(C_theta.*theta) ;
% INHAbs_CLP(end+1)=InhAbs_CL;

end





