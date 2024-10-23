function MBmodel = optimiseMBparams_homeostaticInhibition(MBmodel, PNactivity)
% optimising the "green model"
sigmoid_factor = -1.11;
Sigmoid_derivative = @(x) exp(-x./sigmoid_factor)./((1.+exp(-x./sigmoid_factor)).^2) ;
Sigmoid = @(x) 1./(1.+exp(-x./sigmoid_factor)) ;
 
tic

    % separate optimisation from model creation
    % T=5;
    % theta= abs(normrnd(T,T*(5.6/21.5),[n 1])); %% avoid negative values of theta
    % theta(theta>70)=70;
    % theta(theta<0.01)=0.01;

%n=2000 KCs, numtrainingSamples=15
n = MBmodel.nKCs;
theta = MBmodel.theta;
nResponses = size(PNactivity,2)*size(PNactivity,3);

     %% inhibitory plasticity compensation
            % again, why the multiplication?
        %     thisW_ActivityBasedComp_inhibitionPlast= 5.*thisW;
        % weights = thisW_ActivityBasedComp_inhibitionPlast;%my code
    weights = MBmodel.PNtoKC;
    C_theta = 1;
            % cw=1;
            % Activations =zeros(n,odors*numtrainingSamples);
            % ActivationsDummy= zeros(n, odors*numtrainingSamples);
    APLgains_array = zeros(2000,1);
    A0 = 0.51; %setpoint for KC lifetime sparseness (each KC 
    % should show activity half the time)
    epsilon = A0 * 0.06;       
    eta_APL_1 = 0.02;
    eta_APL_2 = 0.0000001;

        % t=1;
    conditions = false;
    nLoops = 0;
            % A=zeros(n,odors*numtrainingSamples);
            % Y_disInh=zeros(n,odors*numtrainingSamples);
            % codingLevelDummy=[];     
    detectDeadEnd = true;
    if detectDeadEnd
        APLgain_prev = APLgains_array;
        C_prev = C_theta;
    end
    while(~conditions)
        nLoops = nLoops+1;
        if nLoops>1e3
            error('Optimisation did not converge')
        end
                % A = multiprod( ( thisW_ActivityBasedComp_inhibitionPlast)',PNtrials(:,:,1:numtrainingSamples));
                % A=reshape(A,n,odors*numtrainingSamples);
        A = weights' * PNtrials;
                % Y_disInh=(( A-(C_theta.*theta_inhibitionPlast))>0 ).*( A-(C_theta.*theta_inhibitionPlast));
                % codingLevelDummy=  (sum(Y_disInh>0,1)/n);
        Y_disInh = A - C_theta * theta;
        % Y_disInh(Y_disInh<0) = 0;%probably not necessary/right
                % CL_disInh=mean(codingLevelDummy);
        CL_disInh = mean(Y_disInh(:)>0);
        
                % depsi1_dy=(exp(1.*Y_disInh)./((1+exp(1.*Y_disInh)).^2));
                % depsi1_dy(isnan(depsi1_dy))=0;
        depsi_dy = Sigmoid_derivative(Y_disInh);
                % depsi1_dtheta= -(Y_disInh>0).* depsi1_dy.* (repmat(theta_inhibitionPlast,1,odors*numtrainingSamples));
        depsi_dtheta = -(Y_disInh>0) .* depsi_dy .* repmat(theta,1,nResponses);
                % Grad= (((CL_disInh)-0.19)*(1/(n*odors*numtrainingSamples))*(sum(depsi1_dtheta(:) ) ) );
        Grad = (CL_disInh-0.20) * mean(depsi_dtheta(:));
                
        C_theta = C_theta - (mean(Y_disInh(:)) * Grad);
        
        if (C_theta<0 )
            error('the scale factor in green model is -ve')
        end
                
            % Activations = multiprod( ( thisW_ActivityBasedComp_inhibitionPlast)',PNtrials(:,:,1:numtrainingSamples));
            % Activations= reshape(Activations,n,odors*numtrainingSamples);
            
            % Y_incInh= (Activations-(repmat(APLgains_model6,1,odors*numtrainingSamples).*(repmat(sum(Activations),n,1) ))-(C_theta.*theta_inhibitionPlast)>0 ).*( Activations-(repmat(APLgains_model6,1,odors*numtrainingSamples).*(repmat(sum(Activations),n,1)))-(C_theta.*theta_inhibitionPlast));
        Y_incInh = A - repmat(APLgains_array,1,nResponses) .* repmat(sum(A,1),n,1) - C_theta .* theta;
        CL_incInh=mean(Y_incInh(:)>0);

        avgAKcs = mean(Y_incInh,2);
        errorInActivity = avgAKcs-A0;
                % dsig_dy=(exp(1.*Y_incInh)./((1+exp(1.*Y_incInh)).^2));
                % dsig_dy(isnan(dsig_dy))=0;
        dsig_dy = Sigmoid_derivative(Y);
        if any(isnan(dsig_dy(:)))
            fprintf('found NaNs in iteration %d',nLoops)
        end
        dsig_dy(isnan(dsig_dy))=0;
    
                    % xk_=repmat(sum(Activations),n,1);
                    % d_yj_alphaj= (-1.*mean(xk_,2));
        d_yj_alphaj = -1. * mean( sum(A,1) );
                    % dAct_dalpha= sum(Activations,1);
        dAct_dalpha = sum(A,1);
                    % dsig_dalpha= -(Y_incInh>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
                    % dYik_dalphai= -(repmat(dAct_dalpha,n,1));%not used...
        dsig_dalpha = -(Y_incInh>0).*(repmat(dAct_dalpha,n,1)).*  dsig_dy;
        APLgains_array = APLgains_array - eta_APLgain*(Grad_alpha);
                
        % for some KCs change in alpha for activity equalization can cancel out
        % this happen in few KCs towards approaching the global minima for
        % all KCs, the to avoid stucking in a local minima, we add a small noise
        % on top of the calculated gradient.
                
        Grad_alpha1 = eta_APL_1 .* (CL_incInh-0.10) .* mean(dsig_dalpha,2)  + 0.0001*randn(n,1);
        
        Grad_alpha2 = eta_APL_2 * errorInActivity .* d_yj_alphaj + 0.0001*randn(n,1);
                
        Grad_alpha= Grad_alpha1+Grad_alpha2 ;
                
        APLgains_array = APLgains_array - 0.001.*Grad_alpha;
               
        % sanity check
        if detectDeadEnd
            if all(APLgains_array==APLgain_prev) & (C_theta==C_prev)
                fprintf('loop %d',nLoops)
                error('Values remained unchanged without fulfilling the criteria!')
            end
            APLgain_prev = APLgains_array;
            C_prev = C_theta;
        end

                
        %% check constraints
        Y_disInh = A - C_theta * theta;
        CL_disInh = mean(Y_disInh(:)>0);
                    
        Y_incInh = A - repmat(APLgains_array,1,nResponses) .* repmat(sum(A,1),n,1) - C_theta .* theta;
        CL_incInh=mean(Y_incInh(:)>0);
    
        avgAKcs = mean(Y_incInh,2);
                    
        conditions= all(abs(avgAKcs-A0)<epsilon)  & ( abs(round(CL_incInh,3)-0.10) <=0.015 ) & ( abs( CL_disInh/CL_incInh - 2.0) <0.3 );
                
        nLoops = nLoops+1;
    end
            
    disp(sprintf('Optimisation took %d loops',nLoops))

MBmodel.C_theta = C_theta;
MBmodel.alpha = APLgains_array;
MBmodel.CL_disInh = CL_disInh;
MBmodel.CL_incInh = CL_incInh;

end