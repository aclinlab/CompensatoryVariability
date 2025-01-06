%% homeostatic adaptation impacts KC odor responses in numerical simulation

prompt = "Do you to generate MB models automatically? Y/N [N]: ";
txt = input(prompt,"s");
if upper(txt)=='Y'
    ModelGenerationOnly_Philippe; % does the job for now
end

%make these global to simplify use of local function (yeah, bad style...)
global saveFigures
global savepath
global saveFormats
saveFigures = true
saveFormats = {'.fig','.png','.svg'};
savepath = './figures/';

% clean up further at some point, separating model generation from
% olfactory response simulation

% reorder odor trials such that the same odors are next to each other
PNtrials_reordered = reshape(permute(reshape(PNtrials,size(PNtrials,1),[],numTrials),[1 3 2]),size(PNtrials,1),[]);
K_odors = size(PNtrials_reordered,2)/numTrials
% populationResponses = reshape(permute(reshape(populationResponses,size(populationResponses,1),100,[]),[1 3 2]),size(populationResponses,1),[]);

%% calculate the cumulated KC population response to PN input
populationResponses = {};
for ii=1:length(MBmodels_allVar)
    y = calculateKCresponse(MBmodels_allVar{ii}, reshape(PNtrials_reordered,size(PNtrials_reordered,1),[]));
    popResp = sum(y,1);
    populationResponses{ii} = popResp;
end
populationResponses = cell2mat(populationResponses');

meanResp_perMB = nan(length(MBmodels_allVar), size(PNtrials_reordered,2)/numTrials);
for kk=1:size(meanResp_perMB,2)
    meanResp_perMB(:,kk) = mean( populationResponses(:, (kk-1)*numTrials+1 : kk*numTrials ),2);
end

% potentially transform into a GCaMP signal estimate


%% show the results of original model
fig = figure('Position', [10 500 900 200]);
% imagesc(populationResponses(:,1:odors))
imagesc(populationResponses/1000)
xlabel('odor No (repeated trials)')
ylabel('MB model id')
title(sprintf('random MB model:\nKC population activity responses to PN input activity patterns ("odors")'))
colorbar
% add grid lines
hold on
y = [gca().YLim];
for k = 1:numTrials:(K_odors *numTrials)
    x = [k k];
    line(x,y,'Color','w','LineStyle','-');
    %To make sure the grid is visible over all pixel colors, I'll use the 
    % trick of superimposing two line objects
    % plot(x,y,'Color','k','LineStyle',':');
end

% this is a local function
save_figure_perhaps(fig, 'KC-population-responses_odorInput_naiveState_heatmap');

% a boxplot of the averaged responses per model (kinda represents an
% experiment analysis)
fig = figure('Position', [10 10 300 750]);
ax=axes();
% boxplot(ax,meanResp_perMB/1000);
odorSet = [1:109, 140:length(meanResp_perMB)];
xlab = {}
% for kk=10:10:K_odors
for kk=odorSet(mod(odorSet,10)==0)
    xlab = [xlab,repmat({''},1,9),{num2str(kk)}];
end
% boxplot(ax, meanResp_perMB(:, odorSet)/1000, 'Positions',odorSet);
% xlabel('odor No (repeated trials)')
% ylabel('odor response (mean of 30 trials)')
% title(sprintf('random MB model:\nKC population activity responses to PN input ("odors")'))
% xticklabels(xlab)

boxplot(ax, meanResp_perMB(:, odorSet)/1000, 'Positions',odorSet,...
    'orientation','horizontal');
ylabel('odor No (repeated trials)')
xlabel('odor response (mean of 30 trials)')
title(sprintf('random MB model:\nKC population activity responses to PN input ("odors")'))
% yticklabels(xlab)
yticklabels(odorList(odorSet))

ax.Box = 'off';

save_figure_perhaps(fig, 'KC-population-responses_odorInput_naiveState_boxplot');

%% based on visual appearance, select some odors (trials of the same "odor") that gives nice responses
% Isoamyl acetate aka isopentyl acetate is odor #94
% selectOdors = 1171:1200;
% numTrials = size(PNtrials,2)/100
% index of odor isoamyl-acetate (=isopentyl-acetate) can be retrieved via
%    [~,~,odorList] = propagateORNs2PNs();
isoamyl_index = 94; 

selectOdors = ((isoamyl_index-1)*numTrials+1):(isoamyl_index*numTrials);

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
fig = figure('Position', [10 50 1500 200]);
% imagesc(populationResponses(:,1:odors))
imagesc(populationResponses_disInh/1000)
xlabel('odor No (repeated trials)')
ylabel('MB model id')
title(sprintf('dis-inhibited random MB model:\nKC population activity responses to PN input activity patterns ("odors")'))
colorbar
% add grid lines
hold on
y = [gca().YLim];
for k = 1:numTrials:(K_odors*numTrials)
    x = [k k];
    plot(x,y,'Color','w','LineStyle','-');
    %To make sure the grid is visible over all pixel colors, I'll use the 
    % trick of superimposing two line objects
    % plot(x,y,'Color','k','LineStyle',':');
end

save_figure_perhaps(fig,'KC-population-responses_odorInput_disinhibited')

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
fig = figure;
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

save_figure_perhaps(fig, 'KCresp_normalVSdisinhibited_3Cthetas')

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

%% plot this into the same figure

lines_redCthtea = plot([zeros(1,length(selectOdors)); ones(1,length(selectOdors))], ...
   [mean(populationResponses_redCtheta(:,selectOdors),1); mean(populationResponses_redCtheta_disInh(:,selectOdors),1)] ...
   ,'go-' ...
   );
%ultimately, I had to hard-code the legend... MatLab sucks at this -_-
lgd = legend([lines_origCthtea(1), lines_incCthtea(1),lines_redCthtea(1)], 'orig C_{theta}',sprintf('%.2f C_{theta}',round(MBmodels_allVar_adapted{1}.C_theta/MBmodels_allVar{1}.C_theta,1)), sprintf('%.2f C_{theta}',round(MBmodels_allVar_redCtheta{1}.C_theta/MBmodels_allVar{1}.C_theta,1)) );


%% prepare a systematic heatmap
% define the range of values for all variants of this

% alphaFactors = 0.0:0.025:2;
% thetaFactors = 0.0:0.04:2;

% weightsFactors = 0.:0.02:2;
% thetaFactors = 0.0:0.04:2; 

% weightsFactors = 0.0:0.025:2;
% alphaFactors = 0.0:0.04:2;
% thetaFactors = 0.0:0.05:2;

weightsFactors = 0.0:0.01:2;
alphaFactors = 0.0:0.02:2;
thetaFactors = 0.0:0.025:2;

%weightsFactors = 0.2:0.2:6;
%alphaFactors = 0.0:0.05:3;

% weightsFactors_extended = 0.0:0.1:3;
% alphaFactors_extended = 0.0:0.04:6;
% thetaFactors_extended = 0.0:0.1:4;


%% Compare the effect of changing threshold and changing alpha 

testedOdors = PNtrials_reordered(:,selectOdors);
compare_Ctheta_alpha = nan(length(alphaFactors),length(thetaFactors));

for aa=1:length(alphaFactors)  %multiplier for alpha (APL gain)
    for tt=1:length(thetaFactors)  %mutipliers for C_theta (spike threshold)
        MBmodels_allVar_mod = MBmodels_allVar;
        for ii=1:length(MBmodels_allVar_mod)
            %adjust model parameters
            MBmodels_allVar_mod{ii}.C_theta = thetaFactors(tt)*MBmodels_allVar_mod{ii}.C_theta;
            MBmodels_allVar_mod{ii}.alpha = alphaFactors(aa)*MBmodels_allVar_mod{ii}.alpha;

            y = calculateKCresponse(MBmodels_allVar_mod{ii}, testedOdors );
            popResp = sum(y,1);
            compare_Ctheta_alpha(aa,tt) = mean(popResp);
        end

    end
end
%% keep the plotting separate because the previous computation takes a while
fig = figure;
% imagesc(alphaFactors, thetaFactors, singleOdor_popResp_2D)
[X,Y] = meshgrid(alphaFactors, thetaFactors);
% surface(X, Y, compare_Ctheta_alpha', 'edgeColor','none')
contourf(X,Y,compare_Ctheta_alpha')

xlabel('alpha increasing')
ylabel('C_{theta} increasing')
colorbar
title('total KC response for a given odor (mean of 30 trials) depending on APL gain and spike threshold')

save_figure_perhaps(fig, 'APLgain-Ctheta_totalKCresp_1odor_heatmap')

%% explore the responses more systematically:
% -> how much does KC activity change if we take out APL feedback, depending on alpha?
%  -> this corresponds to a horizontal traverse of the previous colormap
%  offsetting by the Y_{alpha=0} activity level

fig = figure;
legentries = {};
% for givenCtheta=[0.5,0.7,1,1.3,1.5]
for givenWeight=[0.75,1,1.3,1.5]
    [~,givenWeightIdx] = min(abs(thetaFactors-givenWeight));
    y_alpha0 = compare_Weights_alpha(givenWeightIdx,1);
    
    plot(alphaFactors, y_alpha0 - compare_Weights_alpha(givenWeightIdx,:))
    hold on
    legentries{end+1} = sprintf('PN->KC strength= %.1f',givenWeight);
end
title('How much does KC activity change by APL>Ort histamine application')
xlabel('alpha (APL feedback strength)')
ylabel('Delta KC activity given \alpha vs \alpha=0')

legend(legentries)

save_figure_perhaps(fig, 'deltaKC_causedByDisinhibition_varWeights');

%% Compare the effect of reducing weights and increasing C_theta 
% presumably the effect is pretty much the same, but with characteristic
% dynamics?

testedOdors = PNtrials_reordered(:,selectOdors);
compare_Weights_Ctheta = nan(length(weightsFactors),length(thetaFactors));

for ww=1:length(weightsFactors)  %multiplier for alpha (APL gain)
    for tt=1:length(thetaFactors)  %mutipliers for C_theta (spike threshold)
        MBmodels_allVar_mod = MBmodels_allVar;
        for ii=1:length(MBmodels_allVar_mod)
            %adjust model parameters
            MBmodels_allVar_mod{ii}.C_theta = thetaFactors(tt)*MBmodels_allVar_mod{ii}.C_theta;
            MBmodels_allVar_mod{ii}.PNtoKC = weightsFactors(ww)*MBmodels_allVar_mod{ii}.PNtoKC;

            y = calculateKCresponse(MBmodels_allVar_mod{ii}, testedOdors );
            popResp = sum(y,1);
            compare_Weights_Ctheta(ww,tt) = mean(popResp);
        end

    end
end

%% plot the resulting 2D-map
fig = figure;
[X,Y] = meshgrid(weightsFactors, thetaFactors);
surface(X,Y, compare_Weights_Ctheta','edgeColor','none')
xlabel('weights increasing')
ylabel('C_{theta} increasing')
colorbar
title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and spike threshold')

if saveFigures
    for format=saveFormats
        saveas(fig,[savepath, 'ExcWeights-Ctheta_totalKCresp_1odor_heatmap', format{:}])
    end
end

% same as a contour plot
fig = figure;
% imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
contourf(X,Y, compare_Weights_Ctheta', 15)
xlabel('weights increasing by factor')
ylabel('C_{\theta} increasing by factor')
colorbar
title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')

save_figure_perhaps(fig, 'ExcWeights-Ctheta_totalKCresp_1odor_contour');

%% Compare the effect of reducing weights and changing alpha 
% presumably the effect is pretty much the same, but with characteristic
% dynamics?
testedOdors = PNtrials_reordered(:,selectOdors);
compare_Weights_alpha = nan(length(weightsFactors),length(alphaFactors), length(MBmodels_allVar_mod));

for ww=1:length(weightsFactors)  %multiplier for alpha (APL gain)
    for aa=1:length(alphaFactors)  %mutipliers for C_theta (spike threshold)
        MBmodels_allVar_mod = MBmodels_allVar;
        for ii=1:length(MBmodels_allVar_mod)
            %adjust model parameters
            MBmodels_allVar_mod{ii}.PNtoKC = weightsFactors(ww)*MBmodels_allVar_mod{ii}.PNtoKC;
            MBmodels_allVar_mod{ii}.alpha = alphaFactors(aa)*MBmodels_allVar_mod{ii}.alpha;

            y = calculateKCresponse(MBmodels_allVar_mod{ii}, testedOdors );
            popResp = sum(y,1);
            compare_Weights_alpha(ww,aa,ii) = mean(popResp);
        end

    end
end
compare_Weights_alpha = mean(compare_Weights_alpha,3);
%% plot the resulting 2D-map

fig = figure;
% imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
[X,Y] = meshgrid(weightsFactors,alphaFactors);
surface(X,Y, compare_Weights_alpha', 'edgeColor','none')
xlabel('weights increasing by factor')
ylabel('alpha increasing by factor')
colorbar
title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')

save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_heatmap')

% same as a contour plot
fig = figure;
% imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
[X,Y] = meshgrid(weightsFactors,alphaFactors);
contourf(X,Y, compare_Weights_alpha'/1000, 25)
xlabel('weights multiplied by factor')
ylabel('alpha multiplied by factor')
colorbar
title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')

save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_contour');

%% calculate this in 3-D: modify weights, alpha and theta

testedOdors = PNtrials_reordered(:,selectOdors);

[w,a,th] = meshgrid( alphaFactors, weightsFactors, thetaFactors);
% [w,a,th] = ndgrid(weightsFactors, alphaFactors, thetaFactors);
threeParamMap_weightAlphaTheta = nan(length(weightsFactors),length(alphaFactors),length(thetaFactors), length(MBmodels_allVar));

MBmodels_allVar_mod = MBmodels_allVar;
for ii=1:length(MBmodels_allVar_mod)
    for ww=1:length(weightsFactors)  %multiplier for alpha (APL gain)
        MBmodels_allVar_mod{ii}.PNtoKC = weightsFactors(ww)*MBmodels_allVar{ii}.PNtoKC;
        for aa=1:length(alphaFactors)  %mutipliers for C_theta (spike threshold)
            MBmodels_allVar_mod{ii}.alpha = alphaFactors(aa)*MBmodels_allVar{ii}.alpha;
            for tt=1:length(thetaFactors)
                MBmodels_allVar_mod{ii}.C_theta = thetaFactors(tt)*MBmodels_allVar{ii}.C_theta;
                y = calculateKCresponse(MBmodels_allVar_mod{ii}, testedOdors );
                popResp = sum(y,1);
                threeParamMap_weightAlphaTheta(ww,aa,tt,ii) = mean(popResp);
            end
        end
    end
end
compareAllParams = mean(threeParamMap_weightAlphaTheta,4);


%% Make a 3D plot of iso-surfaces 

fig = figure;
ax = axes;
for level=(1:2:14)*1e4
    isosurface(w,a,th,compareAllParams, level);
    % isosurface(compareAllParams, level);
end
xlabel('weights increasing by factor')
ylabel('alpha increasing by factor')
zlabel('C_{\theta} increasing')
for ch=[1:length(ax.Children)-2, length(ax.Children)]
    ax.Children(ch).FaceAlpha = 0.4;
end
view(3)

save_figure_perhaps(fig, 'threeParams_totalKCresp_1odor_isosurface');



%% What proportion of KCs is active?

% [w,a,th] = meshgrid( alphaFactors, weightsFactors, thetaFactors);
activeKCproportion_weightAlpha = nan(length(weightsFactors),length(alphaFactors), length(MBmodels_allVar));

MBmodels_allVar_mod = MBmodels_allVar;
for ii=1:length(MBmodels_allVar_mod)
    for ww=1:length(weightsFactors)  %multiplier for alpha (APL gain)
        MBmodels_allVar_mod{ii}.PNtoKC = weightsFactors(ww)*MBmodels_allVar{ii}.PNtoKC;
        for aa=1:length(alphaFactors)  %mutipliers for C_theta (spike threshold)
            MBmodels_allVar_mod{ii}.alpha = alphaFactors(aa)*MBmodels_allVar{ii}.alpha;
    
            y = calculateKCresponse(MBmodels_allVar_mod{ii}, testedOdors );
            activeProp = mean(y>0,1);
            activeKCproportion_weightAlpha(ww,aa,ii) = mean(activeProp);
        end
    end
end

figure
[X,Y] = meshgrid(weightsFactors(1:find(weightsFactors>1.15,1)) , alphaFactors(1:find(alphaFactors>1.15,1)) );
% surface(X,Y, mean(activeKCproportion_weightAlpha, 3)', 'edgeColor','none')
% surface(X,Y, mean(activeKCproportion_weightAlpha(1:find(weightsFactors>1.15,1), 1:find(alphaFactors>1.15,1)), 3)', 'edgeColor','none')
contourf(X,Y, mean(activeKCproportion_weightAlpha(1:find(weightsFactors>1.15,1), 1:find(alphaFactors>1.15,1)), 3)', 20)
xlabel('weights increasing by factor')
ylabel('alpha increasing by factor')
colorbar
title('proportion of active KCs for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')



%% mimic bar plots from Fig 5, 
% i.e. get some odor responses to normal condition and disinhibited

reductionFactors = struct('alphaExcIn',0.9, 'alphaAPLgain',0.25,...
                          'alPrimeExcIn', 0.85, 'alPrimeAPLgain', 0.25,...
                          'gammaExcIn', 0.65, 'gammaAPLgain', 0.25)

barplotResponses = struct();

% unmodified (vanilla) model responses
barplotResponses.vanillaResponses = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',1,'ExcIn',1));
% disinhibit model with original parameters
barplotResponses.vanillaDisinh = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',0,'ExcIn',1));

% modify responses of model by applying excitory-weights and alpha factors
barplotResponses.alphaLobeAdapted = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',reductionFactors.alphaAPLgain, ...
                                   'ExcIn',reductionFactors.alphaExcIn) );
% disinhibit the modified model
barplotResponses.alphaLobeDisinh = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',0, ...
                                   'ExcIn',reductionFactors.alphaExcIn) );

% modify responses of "alPrime lobe" (alpha'/beta') by applying excitory-weights and alpha factors
barplotResponses.alPrimeLobeAdapted = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',reductionFactors.alPrimeAPLgain, ...
                                   'ExcIn',reductionFactors.alPrimeExcIn) );
% disinhibit the modified "alPrime lobe" model
barplotResponses.alPrimeLobeDisinh = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',0, ...
                                   'ExcIn',reductionFactors.alPrimeExcIn) );

% modify responses of "gamma lobe" by applying excitory-weights and alpha factors
barplotResponses.gammaLobeAdapted = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',reductionFactors.gammaAPLgain, ...
                                   'ExcIn',reductionFactors.gammaExcIn) );
% disinhibit the modified "gamma lobe" model
barplotResponses.gammaLobeDisinh = get_modified_MBmodel_responses(MBmodels_allVar,...
                testedOdors,struct('APLgain',0, ...
                                   'ExcIn',reductionFactors.gammaExcIn) );

%% alpha lobe figure
[fig, ax] = make_fig5like_barplot(barplotResponses.vanillaResponses, ...
                     barplotResponses.vanillaDisinh, ...
                     barplotResponses.alphaLobeAdapted, ...
                     barplotResponses.alphaLobeDisinh );
title('my guess for alpha lobe')
save_figure_perhaps(fig, 'barPlot_totalKCresp_mimicFig5_alphaLobe');

%% alPrime lobe figure
[fig, ax] = make_fig5like_barplot(barplotResponses.vanillaResponses, ...
                     barplotResponses.vanillaDisinh, ...
                     barplotResponses.alPrimeLobeAdapted, ...
                     barplotResponses.alPrimeLobeDisinh );
title("my guess for alpha' lobe")
save_figure_perhaps(fig, 'barPlot_totalKCresp_mimicFig5_alPrimeLobe');

%% gamma lobe figure
[fig, ax] = make_fig5like_barplot(barplotResponses.vanillaResponses, ...
                     barplotResponses.vanillaDisinh, ...
                     barplotResponses.gammaLobeAdapted, ...
                     barplotResponses.gammaLobeDisinh );
title('my guess for gamma lobe')
save_figure_perhaps(fig, 'barPlot_totalKCresp_mimicFig5_gammaLobe');

%% local function definitions


function responses = get_modified_MBmodel_responses(MBmodels_allVar, testedOdors, modiFactors)
    %% requires instances of MBmodels, odorInput (24xK matrix), and a struct 
    % of multiplication factors, containing fields "APLgain" and "ExcIn"
    responses = nan(length(MBmodels_allVar),size(testedOdors,2));
    MBmodels_allVar_modif = MBmodels_allVar;
    for ii=1:length(MBmodels_allVar_modif)
        MBmodels_allVar_modif{ii}.alpha = modiFactors.APLgain *...
                            MBmodels_allVar_modif{ii}.alpha;
        MBmodels_allVar_modif{ii}.PNtoKC = modiFactors.ExcIn * ...
                            MBmodels_allVar_modif{ii}.PNtoKC;
        y = calculateKCresponse(MBmodels_allVar_modif{ii}, testedOdors );
        popResp = sum(y,1);
        responses(ii,:) = popResp;
    end
end

function [fig,ax] = make_fig5like_barplot(vanillaResponses, vanillaDisinh, ...
                            adaptedModelResp, adaptedModelDisinh)
    
    disinhibitionColor = [1,0.498, 0.0549];
    normalColor = [0.7, 0.7, 0.7];
    numModels = size(vanillaResponses,1);
    
    fig = figure;
    ax = axes;
    barData = mean( [ mean(vanillaResponses,2) , ...
        mean(vanillaDisinh,2) , ...
        mean(adaptedModelResp,2) , ...
        mean(adaptedModelDisinh,2) ] ) / 1000;
    b = bar(reshape(barData,2,2)');
    b(2).FaceColor = disinhibitionColor ;
    b(1).FaceColor = normalColor;
    ylabel('total KC activity [A.U.]', 'fontsize',16)
    xticklabels({'naive','adapted'})
    ax.XAxis.FontSize=16;
    
    hold on
    plot(repmat([0.85;1.15],1,numModels), [mean(vanillaResponses,2) , ...
        mean(vanillaDisinh,2)]'/1000 , 'color','#606060')
    plot(repmat([1.85;2.15],1,numModels), [mean(adaptedModelResp,2) , ...
        mean(adaptedModelDisinh,2)]'/1000 , 'color','#606060')
    
    legend({'normal','disinhibited'}, 'fontsize',16)
end

function save_figure_perhaps(fig, figureName)
    global saveFigures
    global savepath
    global saveFormats
    if saveFigures
        for format=saveFormats
            %use switch with print instead of savepath for finer control
            switch format{:}
                case '.png'
                    print([savepath, figureName, format{:}], '-dpng')
                case '.jpg'
                    print([savepath, figureName, format{:}], '-djpeg')
                case '.svg'
                    print('-vector', [savepath, figureName, format{:}], '-dsvg')
                case '.fig'
                    saveas(fig,[savepath, figureName, format{:}])
            end
        end
    end
end