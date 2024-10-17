function MBmodel = buildMBmodel(nKCs, nPNs ,isVar_Nclaws, isVar_weights, ...
                                    isVar_Theta, varargin)

p = inputParser;
addRequired(p,'nKCs', @isnumeric);
addRequired(p,'nPNs', @isnumeric);

addRequired(p,'isVar_Nclaws');
addParameter(p,'Nclaw_mean', 6, @isnumeric);
addParameter(p,'Nclaw_std', 1.7, @isnumeric);
addParameter(p,'Nclaw_limits', [2,11], @isnumeric);

addRequired(p,'isVar_weights'); % PN to KC weights

addRequired(p,'isVar_Theta'); %about spiking thresholds
addParameter(p,'theta_mean', 10, @isnumeric);
addParameter(p,'theta_std', 10*5.6/21.5, @isnumeric);
addParameter(p,'theta_limits', [0.01,70], @isnumeric);

% parse(p, 'nKCs', 'nPNs', 'isVar_Nclaws', 'isVar_weights', 'isVar_Theta', varargin);
parse(p, nKCs, nPNs, isVar_Nclaws, isVar_weights, isVar_Theta, varargin{:});

% start making the MBmodel struct
MBmodel = p.Results;

% make the number of claws, distinguish homogeneous vs variable case
if p.Results.isVar_Nclaws
    %select number of claws randomly 
    % -> for models depicted with green in Fig2
    clawsNo=(normrnd(p.Results.theta_mean, p.Results.theta_std, [1, MBmodel.nKCs])); 
    clawsNo(clawsNo < p.Results.Nclaw_limits(1)) = p.Results.Nclaw_limits(1);
    clawsNo(clawsNo > p.Results.Nclaw_limits(2)) = p.Results.Nclaw_limits(2);
    clawsNo = round(clawsNo);
else
    clawsNo = MBmodel.theta_mean*ones(MBmodel.nKCs,1);
end
% set up individual KC connectivity, sample either variable or
% homogenously in terms of PN-KC connections
for ii=1:MBmodel.nKCs %iterate through KCs
    %draws correspond to PN ids
    PnToKc{ii} = randsample(MBmodel.nPNs, clawsNo(ii), true); 
end

%initialize the weights matrix between KC and PN
weights = zeros(MBmodel.nPNs, MBmodel.nKCs);
if p.Results.isVar_weights
    for jj=1:MBmodel.nKCs
        % for ii=1:length(PnToKc{k})
        upstreamPNs = PnToKc{jj};
        thisWeight = exp(-0.0507+0.3527*randn(length(upstreamPNs),1)); %randn=Norm(0,1, size)
        weights(upstreamPNs, jj) = thisWeight;
    end
else
    for jj=1:MBmodel.nKCs
        upstreamPNs = PnToKc{jj};
        weights(upstreamPNs, jj) = 1;
    end
end

if p.Results.isVar_Theta
    theta = p.Results.theta_mean + p.Results.theta_std*randn(MBmodel.nKCs,1);
    theta(theta<p.Results.theta_limits(1)) = p.Results.theta_limits(1);
    theta(theta>p.Results.theta_limits(2)) = p.Results.theta_limits(2);
else
    theta = p.Results.theta_mean*ones(MBmodel.nKCs,1);
end

MBmodel.theta = theta;
MBmodel.PNtoKC = weights;


        