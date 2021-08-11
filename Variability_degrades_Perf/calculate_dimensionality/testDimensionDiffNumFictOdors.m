%% This script reproduces Fig S2A - calculating dimensionality with different numbers of fictitious odors

n = 2000;
numFlies = 50;
numFictOdors_values = [10 20 50 100 200 500 1000 2000 5000 10000 20000 50000];% 100000];
numTrials = 10;
dim = zeros(length(numFictOdors_values), numTrials);

load('hallem_olsen.mat');

PN = hallem_olsen(1:110,:)';
% get 1000 odors from the 110 PNs
for Pn=1:24
    [prob,bins]=hist(PN(Pn,:),100);
    prob=prob/sum(prob);
    binsAll(Pn,:) = bins;
end
%% recover the rescaling factors from the original PNs used to tune the coding level
% replace this line with the .mat file storing the variable 'PNtrials'
% which provided the inputs used to tune the coding level of the model for
% which you are measuring the dimensionality
load('CL_performance_APLnotzero_random_homog100_0.1_1.mat'); 
PNs_1_=PNtrials(:,:,1);
% get the maximum bin center for each PN derived from the original
% Hallem-Olsen data
maxRespBinPerPN = max(binsAll,[],2);

% get the maximum response in the rescaled randomly resampled PNs
maxRespPNsRescaledPerPN = max(PNs_1_,[],2);

PNsAboveBestFit = true(24,1);

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

while sum(PNsAboveBestFit)
    sum(PNsAboveBestFit)
    
    p = polyfit(maxRespBinPerPN(PNsAboveBestFit), maxRespPNsRescaledPerPN(PNsAboveBestFit),1);
    
    % the correct PNs will be above the best fit line because the incorrect PNs
    % are outliers that drag down the line of best fit
    PNsAboveBestFit = maxRespPNsRescaledPerPN > ( maxRespBinPerPN*p(1) + p(2) +0.000001);
    % the 0.000001 is for rounding errors, otherwise you end up in endless
    % loops
    
end

%% plot the data and the best fit line to confirm it passes through the correct PNs
figure,scatter(maxRespBinPerPN,maxRespPNsRescaledPerPN);
hold on
xrange = [min(binsAll(:)) max(binsAll(:))];
plot(xrange, p(1)*xrange + p(2))
for t = 1:numTrials
    t
    for i = 1:length(numFictOdors_values)
        %% create artificial odors, mulOd odors
        
        mulOd = numFictOdors_values(i);
        
        PNs=zeros(24,mulOd);
        for Pn=1:24
            
            PNs(Pn,:)=randsample(bins,mulOd,'true',prob);
        end
        
        % add noise
        PNs = PNs + getPNStdevBhandawat(PNs).*randn(size(PNs));
        
        
        
        
        %% apply the original scaling to the new fictitious odors
        PNs = PNs*p(1) + p(2);
        PNs(PNs<0) = 0;
        
        %% simulate activity for fictitious odors, calculate dimensionality
        numFictOdors_values(i)
        tic
        filename = strcat('CL_performance_APLnotzero_random_homog100_0.1_1.mat');
        load(filename);
        
        homog_resp = getKCactivity(thisW_HomogModel, PNs, APLgainP(2), repmat(thetaH_Ftheta,n,1));
        disp(nnz(homog_resp)/(size(PNs,2)*size(PNs,3)*n));
        
        toc
        
        tic
        
        dim(i,t) = dimInputCurrent(Calc_C(homog_resp));
        toc
    end
end
save('dim_diffNumFictOdors_wAPL_CL0.1.mat','dim','numFictOdors_values');
figure
errorbar(numFictOdors_values, mean(dim,2), ConfInt95(dim,2),'-ok','MarkerFaceColor','k','MarkerSize',4)
set(gca,'units','centimeters','position',[4 4 4 4])
set(gca,'XScale','log')
xlabel('# simulated odors')
ylabel('Calculated dimensionality')
set(gca,'FontSize',12)
set(gca,'LabelFontSizeMultiplier',1)
set(gca,'TickLength',[.03 .03],'TickDir','out')
set(gca,'XTick',10.^(1:5))
% ylim([0.5 0.75])
xlim(10.^[0.5 5])
box off


function Y = getKCactivity(w, PNtrials, APLgain, theta)
odors = size(PNtrials,2);
numTrials = size(PNtrials,3);
n = size(w,2);
[Activations,Y]=deal(zeros(n,odors*numTrials));
for trial = 1:(odors*numTrials)
    
    Activations(:,trial) = w'*PNtrials(:,trial );
    %     resp = Activations(:,trial)-APLgain*repmat(sum(Activations(:,trial),1),n,1)-theta;
    %     resp(resp<0) = 0;
    %     Y(:,trial)=resp;
    Y(:,trial)=(( Activations(:,trial)-APLgain*repmat(sum(Activations(:,trial),1),n,1)-theta)>0 ).*( Activations(:,trial)-APLgain*repmat(sum(Activations(:,trial),1),n,1)-theta);
end

end
