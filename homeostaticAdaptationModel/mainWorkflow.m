%% homeostatic adaptation impacts KC odor responses in numerical simulation

prompt = "Do you to generate MB models automatically? Y/N [N]: ";
txt = input(prompt,"s");
if upper(txt)=='Y'
    ModelGenerationOnly_Philippe; % does the job for now
end

% clean up further at some point, separating model generation from
% olfactory response simulation

% reorder odor trialssuch that the same odors are next to each other
PNtrials_reordered = reshape(permute(reshape(PNtrials,size(PNtrials,1),100,[]),[1 3 2]),size(PNtrials,1),[]);
% populationResponses = reshape(permute(reshape(populationResponses,size(populationResponses,1),100,[]),[1 3 2]),size(populationResponses,1),[]);

%% calculate the cumulated KC population response to PN input
populationResponses = {};
for ii=1:length(MBmodels_allVar)
    y = calculateKCresponse(MBmodels_allVar{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses{ii} = popResp;
end
populationResponses = cell2mat(populationResponses');

% potentially transform into a GCaMP signal estimate

%% show the results of original model
figure('Position', [10 600 900 300]), 
% imagesc(populationResponses(:,1:odors))
imagesc(populationResponses)
xlabel('odor No (repeated trials)')
ylabel('MB model id')
title(sprintf('random MB model:\nKC population activity responses to PN input activity patterns ("odors")'))
colorbar
% add grid lines
hold on
y = [gca().YLim];
for k = 1:numTrials:(odors*numTrials)
    x = [k k];
    plot(x,y,'Color','w','LineStyle','-');
    %To make sure the grid is visible over all pixel colors, I'll use the 
    % trick of superimposing two line objects
    % plot(x,y,'Color','k','LineStyle',':');
end

%% based on visual appearance, select some odors (trials of the same "odor") that gives nice responses
selectOdors = 1171:1200;

singleOdorResp = [mean(populationResponses(:,selectOdors),'all'), 
                  std(populationResponses(:,selectOdors),[],"all") ]; 

%% make a version of these models where there is no APL inhibition (or greatly reduced)
MBmodels_allVar_disInh = MBmodels_allVar;
for ii=1:length(MBmodels_allVar_disInh)
    MBmodels_allVar_disInh{ii}.alpha = 0;
end

populationResponses_disInh = {};
for ii=1:length(MBmodels_allVar_disInh)
    y = calculateKCresponse(MBmodels_allVar_disInh{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses_disInh{ii} = popResp;
end
populationResponses_disInh = cell2mat(populationResponses_disInh');

singleOdorResp_disInh = [mean(populationResponses_disInh(:,selectOdors),'all'), 
                         std(populationResponses_disInh(:,selectOdors),[],'all') ]; 

%% show the results of dis-inhibited original model
figure('Position', [10 50 1500 200])
% imagesc(populationResponses(:,1:odors))
imagesc(populationResponses_disInh)
xlabel('odor No (repeated trials)')
ylabel('MB model id')
title(sprintf('dis-inhibited random MB model:\nKC population activity responses to PN input activity patterns ("odors")'))
colorbar
% add grid lines
hold on
y = [gca().YLim];
for k = 1:numTrials:(odors*numTrials)
    x = [k k];
    plot(x,y,'Color','w','LineStyle','-');
    %To make sure the grid is visible over all pixel colors, I'll use the 
    % trick of superimposing two line objects
    % plot(x,y,'Color','k','LineStyle',':');
end

%% now arbitrarily INCREASE C_theta to COARSELY mimic the phenotype of KC>TRPA1 activation
% calculate the responses, and make a sequential plot for some of the best
% odors

MBmodels_allVar_adapted = MBmodels_allVar;
for ii=1:length(MBmodels_allVar_adapted)
    MBmodels_allVar_adapted{ii}.C_theta = 2*MBmodels_allVar_adapted{ii}.C_theta;
end

populationResponses_adapted = {};
for ii=1:length(MBmodels_allVar_adapted)
    y = calculateKCresponse(MBmodels_allVar_adapted{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses_adapted{ii} = popResp;
end
populationResponses_adapted = cell2mat(populationResponses_adapted');

% singleOdorResp_disInh = [mean(populationResponses_adapted(:,selectOdors),'all'), 
%                          std(populationResponses_adapted(:,selectOdors),[],'all') ]; 


%% analogous to the above, set APL gain to 0 to mimic APL>Ort+histamine
MBmodels_allVar_adapted_disInh = MBmodels_allVar_adapted;
for ii=1:length(MBmodels_allVar_adapted_disInh)
    MBmodels_allVar_adapted_disInh{ii}.alpha = 0;
end

populationResponses_adapted_disInh = {};
for ii=1:length(MBmodels_allVar_adapted_disInh)
    y = calculateKCresponse(MBmodels_allVar_adapted_disInh{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses_adapted_disInh{ii} = popResp;
end
populationResponses_adapted_disInh = cell2mat(populationResponses_adapted_disInh');

% singleOdorResp_disInh = [mean(populationResponses_adapted_disInh(:,selectOdors),'all'), 
%                          std(populationResponses_adapted_disInh(:,selectOdors),[],'all') ]; 

%% plot the slopes
figure
lines_origCthtea = plot([zeros(1,length(selectOdors)); ones(1,length(selectOdors))], ...
   [mean(populationResponses(:,selectOdors),1); mean(populationResponses_disInh(:,selectOdors),1)] ...
   , 'bo-' ...
   ,'DisplayName', 'orig C_{theta}' ... %'orig C_{theta}'...
   );
hold on
lines_incCthtea = plot([zeros(1,length(selectOdors)); ones(1,length(selectOdors))], ...
   [ mean(populationResponses_adapted(:,selectOdors),1); mean(populationResponses_adapted_disInh(:,selectOdors),1)] ...
   ,'ro-' ...
    ,'DisplayName',sprintf('%.2f C_{theta}',round(MBmodels_allVar_adapted{1}.C_theta/MBmodels_allVar{1}.C_theta,1)) ...
  );

ylabel('population activity to given odor averaged across flies')

xticks([0,1])
xticklabels({sprintf('normal / adapted state'),'APL silenced'})
title('KC population responses to trials of the same odor')
xlim([-0.25,1.25])
lgd = legend([lines_origCthtea(1), lines_incCthtea(1)], 'orig C_{theta}',sprintf('%.2f C_{theta}',round(MBmodels_allVar_adapted{1}.C_theta/MBmodels_allVar{1}.C_theta,1) ) )
lgd.AutoUpdate = 'off';

%% now arbitrarily REDUCE C_theta to COARSELY mimic the phenotype of KC>TRPA1 activation
% calculate the responses, and make a sequential plot for some of the best
% odors

MBmodels_allVar_redCtheta = MBmodels_allVar;
for ii=1:length(MBmodels_allVar_redCtheta)
    MBmodels_allVar_redCtheta{ii}.C_theta = 0.5*MBmodels_allVar_redCtheta{ii}.C_theta;
end

populationResponses_redCtheta = {};
for ii=1:length(MBmodels_allVar_redCtheta)
    y = calculateKCresponse(MBmodels_allVar_redCtheta{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses_redCtheta{ii} = popResp;
end
populationResponses_redCtheta = cell2mat(populationResponses_redCtheta');

% singleOdorResp_disInh = [mean(populationResponses_redCtheta(:,selectOdors),'all'), 
%                          std(populationResponses_redCtheta(:,selectOdors),[],'all') ]; 

% also set APL gain to 0
MBmodels_allVar_redCtheta_disInh = MBmodels_allVar_redCtheta;
for ii=1:length(MBmodels_allVar_redCtheta_disInh)
    MBmodels_allVar_redCtheta_disInh{ii}.alpha = 0;
end

populationResponses_redCtheta_disInh = {};
for ii=1:length(MBmodels_allVar_redCtheta_disInh)
    y = calculateKCresponse(MBmodels_allVar_redCtheta_disInh{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses_redCtheta_disInh{ii} = popResp;
end
populationResponses_redCtheta_disInh = cell2mat(populationResponses_redCtheta_disInh');

% singleOdorResp_disInh = [mean(populationResponses_redCtheta_disInh(:,selectOdors),'all'), 
%                          std(populationResponses_redCtheta_disInh(:,selectOdors),[],'all') ]; 

%% plot this into the same fiure

lines_redCthtea = plot([zeros(1,length(selectOdors)); ones(1,length(selectOdors))], ...
   [mean(populationResponses_redCtheta(:,selectOdors),1); mean(populationResponses_redCtheta_disInh(:,selectOdors),1)] ...
   ,'go-' ...
   );
%ultimately, I had to hard-code the legend... MatLab sucks at this -_-
lgd = legend([lines_origCthtea(1), lines_incCthtea(1),lines_redCthtea(1)], 'orig C_{theta}',sprintf('%.2f C_{theta}',round(MBmodels_allVar_adapted{1}.C_theta/MBmodels_allVar{1}.C_theta,1)), sprintf('%.2f C_{theta}',round(MBmodels_allVar_redCtheta{1}.C_theta/MBmodels_allVar{1}.C_theta,1)) );



%% prepare a systematic heatmap
% -> computation takes a while

alphaFactors = 0.05:0.01:2;
thetaFactors = 0.05:0.01:2;
testedOdors = PNtrials_reordered(:,selectOdors);
singleOdor_popResp_2D = nan(length(alphaFactors),length(thetaFactors));

for aa=1:length(alphaFactors)  %multiplier for alpha (APL gain)
    for tt=1:length(thetaFactors)  %mutipliers for C_theta (spike threshold)
        MBmodels_allVar_mod = MBmodels_allVar;
        for ii=1:length(MBmodels_allVar_mod)
            %adjust model parameters
            MBmodels_allVar_mod{ii}.C_theta = thetaFactors(tt)*MBmodels_allVar_mod{ii}.C_theta;
            MBmodels_allVar_mod{ii}.alpha = alphaFactors(aa)*MBmodels_allVar_mod{ii}.alpha;

            y = calculateKCresponse(MBmodels_allVar_mod{ii}, testedOdors );
            popResp = sum(y,1);
            singleOdor_popResp_2D(aa,tt) = mean(popResp);
        end

    end
end
%% keep the plotting separate because the previous computation takes a while
figure
imagesc(singleOdor_popResp_2D)
xlabel('alpha increasing')
ylabel('C_{theta} increasing')
colorbar
title('total KC response for a given odor (mean of 30 trials) depending on APL gain and spike threshold')

