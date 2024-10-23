function PNactivity = generate_PN_odor_responses()

%candidates for arguments to be integrated into function call
numtrainingSamples=15;
% noiseLevels=[0.5,1,2];
noiseLevel = 1;
K_odors = 200;
numTrials = 30;

	% the odors we used in Fig2, 3 and 4, stored in this .mat file.
	load('Data_submitted_fly_wNoise11.mat','PNtrials');
	% the scaled PNs responses, the base 100 fictitious odors.
	PNs_1_ = PNtrials(:,:,1);
nPNs = size(PNs_1_,1);


if K_odors>100
    load('./hallem_olsen.mat');
	PNresponses_hallem = hallem_olsen(1:110,:)';
	% nPNs = size(PNresponses_hallem,1);
	k=0;
	binsAll=zeros(nPNs,100);

	%% resample Hallem-Olsen dataset to generate new realistic odor responses
	resampledResponses = zeros(nPNs,K_odors-100); %we want to get to K total odors, and will retrieve 100 from another file
	for pn=1:nPNs
		[prob,bins] = hist(PNresponses_hallem(pn,:),100);
		prob = prob/sum(prob);
		
		binsAll(pn,:) = bins;
		resampledResponses(pn,k+1:k+K_odors-100)=randsample(bins,(K_odors-100),'true',prob);
	end

	%% recover the rescaling factors
	% get the maximum bin center for each PN derived from the original
	% Hallem-Olsen data
	maxRespBinPerPN = max(binsAll,[],2);
	% get the maximum response in the rescaled randomly resampled PNs
	maxRespPNsRescaledPerPN = max(PNs_1_,[],2);

	% In this while loop:
	% draw a best fit line comparing the maximum response in the rescaled
	% randomly resample PNs to the maximum bin center for each PN in the
	% original H-O data. Most PNs will match, but in some PNs, by random chance
	% they will not have sampled the top response in 100 odors. These PNs will
	% be below the best fit line, while the "matching" PNs will be above. On
	% the next iteration of the while loop, redraw the best fit using only the
	% PNs that lie above the best fit line from the current iteration
	% continue until none of the PNs are above the best fit line (because they
	% lie on it almost exactly)
	PNsAboveBestFit = true(nPNs,1);
	while sum(PNsAboveBestFit)
		% sum(PNsAboveBestFit)
		p = polyfit(maxRespBinPerPN(PNsAboveBestFit), maxRespPNsRescaledPerPN(PNsAboveBestFit),1);
		
		% the correct PNs will be above the best fit line because the incorrect PNs
		% are outliers that drag down the line of best fit
		PNsAboveBestFit = maxRespPNsRescaledPerPN > ( maxRespBinPerPN*p(1) + p(2) +0.000001);
		% the 0.000001 is for rounding errors, otherwise you end up in endless
		% loops
	end
	PNs_1 = (PNs_1_ - p(2))/p(1);
	% PNs_2=  (resampledResponses*p(1))+p(2);
    % replicate the way it was in Nada's code
	PNs = [PNs_1 resampledResponses]; % now this is the 24 by 200 matrix (24 ORNs 200 odors)
else
	PNs = PNs_1_;
end

PNactivity = zeros(nPNs, K_odors, numTrials);
x = PNs(:,1:K_odors);
PNactivity(:,:,1) = x; % int his setting, it's just all the columns
% get hard-coded glomeruli-specific std-value from Bhandawat2007
noiseStdPerPN = getPNStdevBhandawat(x);
% add some noise to the fictitious odor response
PNactivity(:,:,2:end) = x + ... 
			noiseStdPerPN .* ... 
			noiseLevel .* ...
			randn(nPNs, K_odors, numTrials-1);
PNactivity(PNactivity<0)=0;
%rescaling the PN responses in range from 0 to 5 -> not sure why
PNactivity = (PNactivity - min(PNactivity(:))+0 )/(max(PNactivity(:))-min(PNactivity(:))) * 5;
PNactivity = reshape(PNactivity, size(PNactivity,1),[]); %make 2-D


