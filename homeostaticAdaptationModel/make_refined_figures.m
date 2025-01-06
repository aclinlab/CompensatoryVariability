%% making refined figures only

% load('allThreeParamsMap_weightAlphaTheta__NONhomeostatModels.mat')
% load('allThreeParamsMap_weightAlphaTheta__homeostaticModels.mat')

global saveFigures
global savepath
global saveFormats
saveFigures = false;
saveFormats = {'.fig','.png','.svg'};
savepath = './figures/';

% ADJUST HERE!
MBmodels_allVar = MBmodels_allVar_withHomeo;

PNtrials = allTestedOdors{1};

isoamyl_index = 94; 
selectOdors = ((isoamyl_index-1)*numTrials+1):(isoamyl_index*numTrials);

PNtrials_reordered = reshape(permute(reshape(PNtrials,size(PNtrials,1),[],numTrials),[1 3 2]),size(PNtrials,1),[]);
testedOdors = PNtrials_reordered(:,selectOdors);

% verify
prompt = "Do you to generate MB models automatically? Y/N [N]: ";
disp("You are currently operating on the model type ")
disp(MBmodels_allVar{1}.modelType)
txt = input("Are you sure that you selected the right kind of model and saved data? Y/N [N]: ","s");
if ~(upper(txt)=='Y')
    return
end

% get the slice where C_theta is unchanged, varying only APL and excitation
compare_Weights_alpha = compareAllParams(:,:,find(thetaFactors==1) );

% zoom in on the area below 1 (and take some surroundings)
APL_cutoff = find(alphaFactors<=1.1, 1, "last");
Excit_cutoff = find(weightsFactors<=1.10, 1, "last");
theta_cutoff= find(thetaFactors>=0., 1, "first");
reducedSet = struct('compare_Weights_alpha', compare_Weights_alpha(1:Excit_cutoff, 1:APL_cutoff), ...
                    'alphaFactors',alphaFactors(1:APL_cutoff), ...
                    'weightsFactors',weightsFactors(1:Excit_cutoff), ...
                    'thetaFactors', thetaFactors(theta_cutoff:end-10),...
                    'compare_Weights_theta', compareAllParams(1:Excit_cutoff,find(alphaFactors==1),theta_cutoff:end-10) );

%% Compare the effect of reducing weights and changing alpha 

% fig = figure;
% % imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
% [X,Y] = meshgrid(reducedSet.weightsFactors,reducedSet.alphaFactors);
% surface(X,Y, reducedSet.compare_Weights_alpha', 'edgeColor','none')
% xlabel('weights increasing by factor')
% ylabel('alpha increasing by factor')
% colorbar
% title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')
% 
% save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_heatmap')

% same as a contour plot
fig = figure;
% imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
[X,Y] = meshgrid(reducedSet.weightsFactors,reducedSet.alphaFactors);
% mylevels = [0.0:0.2:4.4]; % fits the non-homeostatic models
mylevels = [0.0:0.5:13]; % fits the homeostatic models
% mylevels = [0.0002:0.4:17];
% mylevels = 10.^([-6:0.2:1.2])
contourf(X,Y, reducedSet.compare_Weights_alpha'/1000, mylevels, 'LineWidth',0.5)
cmap = parula(length(mylevels));
% cmap(1,:) = [.3,.3,.3];
% colormap(cmap)
colorbar
hold on
smallerLevels = [0:0.1:1];
% [m,c] = contourf(X,Y, reducedSet.compare_Weights_alpha'/1000, smallerLevels,...
    % 'w--','FaceAlpha',0);
% contourf(X,Y, reducedSet.compare_Weights_alpha'/1000, 25, "ShowText",true,"LabelFormat","%0.2f")

xlabel('weights multiplied by factor')
ylabel('alpha multiplied by factor')
title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')

save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_contour_reducedSize');


%% one up this by adding line plots

fig = figure;
set(fig,'units','normalized','position',[0.1,0.1,0.8,0.8]);
axAPLconst = subplot(5,5, 1:4);
axExcitconst = subplot(5,5, [10:5:25]);
axMain = subplot(5,5, [6:9,11:14,16:19,21:24]);
hold(axExcitconst,'on')
hold(axAPLconst,'on')

mylevels = [0.0:0.5:13]; % fits the homeostatic models

[X,Y] = meshgrid(reducedSet.weightsFactors,reducedSet.alphaFactors);
[m,c] = contourf(axMain, X,Y, reducedSet.compare_Weights_alpha'/1000, mylevels);
% cmap = parula(length(mylevels));
% colormap(cmap)
cbar = colorbar;

% adding line plots
for transectWeightAt=[0.6,0.7,0.8,0.9,1.0]
    plot(axExcitconst, reducedSet.alphaFactors, ...
         reducedSet.compare_Weights_alpha(find(isEqualFloat(weightsFactors,transectWeightAt)),:)/1000,...
         'DisplayName',num2str(transectWeightAt));
    line(axMain, [transectWeightAt,transectWeightAt], reducedSet.alphaFactors([1,end]), linestyle='--', color='k', linewidth=2)
end

for transectAlphaAt=0.4:0.2:1.0
    plot(axAPLconst, reducedSet.weightsFactors, ...
             reducedSet.compare_Weights_alpha(:,find(isEqualFloat(alphaFactors,transectAlphaAt)))/1000,...
             'DisplayName',num2str(transectAlphaAt));
    line(axMain, reducedSet.weightsFactors([1,end]), [transectAlphaAt,transectAlphaAt],...
        linestyle='--', color='k',linewidth=2)
end
view(axExcitconst, 270,90); % Rotate the right subplot
xlabel(axExcitconst, 'APL gain multiplier')
xlabel(axAPLconst, 'excitation multiplier')
xlabel(axMain, 'excitation multiplier')
ylabel(axMain, 'APL gain multiplier')
ylabel(axExcitconst, 'total KC activity [A.U.]')
ylabel(axAPLconst, 'total KC activity [A.U.]')
%mark transect points lines in contourplot
ylabel(cbar,'total KC activity [A.U]') %label the colorbar ("color axis")
axAPLconst.Position(3) = axMain.Position(3); % make both the same size (nec. bc of colorbar)
fontsize(12,'points')
xlim(axAPLconst,[0,reducedSet.weightsFactors(end)])
xlim(axExcitconst,[0,reducedSet.alphaFactors(end)])

save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_contour_reducedSize_fancyVersion');

%% same, but different layout

fig = figure;
set(fig,'units','normalized','position',[0.1,0.1,0.8,0.8]);
axAPLconst = subplot(2,3, 3);
axExcitconst = subplot(2,3, 6);
axMain = subplot(2,3, [1,2,4,5]);
hold(axAPLconst,'on')
hold(axExcitconst,'on')

mylevels = [0.0:0.5:13]; % fits the homeostatic models

[m,c] = contourf(axMain, X,Y, reducedSet.compare_Weights_alpha'/1000, mylevels);
% cmap = parula(length(mylevels));
% colormap(cmap)
cbar = colorbar;

% weightTransects = [0.6,0.7,0.8,0.9,1.0];
weightTransects = 0.625:0.075:1;
linecolors = copper(length(weightTransects));
linecolors= linecolors(end:-1:1,:);
for transectWeightAt=weightTransects
    plot(axExcitconst, reducedSet.alphaFactors, ...
         reducedSet.compare_Weights_alpha(find(isEqualFloat(weightsFactors,transectWeightAt)),:)/1000,...
         'DisplayName',num2str(transectWeightAt), 'color',linecolors(end,:));
    line(axMain, [transectWeightAt,transectWeightAt], reducedSet.alphaFactors([1,end]), ...
        linestyle='--', color=linecolors(end,:), linewidth=1.5)
    linecolors(end,:) = []; %delete last, this method avoids counter var
end

% alphaTransects = 0.28:0.12:1.0;
alphaTransects = 0.20:0.16:1.0;
linecolors = bone(length(alphaTransects)+2);%truncate extras, avoid white on white
linecolors = linecolors(1:end-2,:);
for transectAlphaAt=alphaTransects
    disp(transectAlphaAt)
    plot(axAPLconst, reducedSet.weightsFactors, ...
             reducedSet.compare_Weights_alpha(:,find(isEqualFloat(alphaFactors,transectAlphaAt)))/1000,...
             'DisplayName',num2str(transectAlphaAt), 'color',linecolors(end,:));
    line(axMain, reducedSet.weightsFactors([1,end]), [transectAlphaAt,transectAlphaAt],...
        linestyle='--', color=linecolors(end,:),linewidth=1.5)
    linecolors(end,:) = []; %delete last, this method avoids counter var
end
% view(axExcitconst, 270,90); % Rotate the right subplot
xlabel(axExcitconst, 'APL gain multiplier')
xlabel(axAPLconst, 'excitation multiplier')
xlabel(axMain, 'excitation multiplier')
ylabel(axMain, 'APL gain multiplier')
ylabel(axExcitconst, 'total KC activity [A.U.]')
ylabel(axAPLconst, 'total KC activity [A.U.]')
%mark transect points lines in contourplot
ylabel(cbar,'total KC activity [A.U]') %label the colorbar ("color axis")
% axAPLconst.Position(3) = axMain.Position(3); % make both the same size (nec. bc of colorbar)
% fontsize(12,'points')
xlim(axAPLconst,[0,reducedSet.weightsFactors(end)])
xlim(axExcitconst,[0,reducedSet.alphaFactors(end)])

save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_contour_reducedSize_fancyVersion2');


%% Compare the effect of changing weights and changing theta

% fig = figure;
% % imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
% [X,Y] = meshgrid(reducedSet.weightsFactors,reducedSet.alphaFactors);
% surface(X,Y, reducedSet.compare_Weights_alpha', 'edgeColor','none')
% xlabel('weights increasing by factor')
% ylabel('alpha increasing by factor')
% colorbar
% title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and APL gain')
% 
% save_figure_perhaps(fig, 'ExcWeights-APLgain_totalKCresp_1odor_heatmap')

% same as a contour plot
fig = figure;
% imagesc(weightsFactors, alphaFactors, compare_Weights_alpha')
[X,Y] = meshgrid(reducedSet.weightsFactors,reducedSet.thetaFactors);
% mylevels = [0.0:0.2:4.4]; % fits the non-homeostatic models
% mylevels = [0.0:0.5:13]; % fits the homeostatic models
% mylevels = [0.0002:0.4:17];
% mylevels = 10.^([-6:0.2:1.2])
mylevels = 20;
contourf(X,Y, squeeze(reducedSet.compare_Weights_theta)'/1000, mylevels, 'LineWidth',0.5)
% cmap = parula(length(mylevels));
% cmap(1,:) = [.3,.3,.3];
% colormap(cmap)
colorbar
hold on
% smallerLevels = [0:0.1:1];
% [m,c] = contourf(X,Y, reducedSet.compare_Weights_alpha'/1000, smallerLevels,...
    % 'w--','FaceAlpha',0);
% contourf(X,Y, reducedSet.compare_Weights_alpha'/1000, 25, "ShowText",true,"LabelFormat","%0.2f")

xlabel('weights multiplied by factor')
ylabel('thresholds multiplied by factor')
title('total KC response for a given odor (mean of 30 trials) depending on excitatory strength and spike threshold')

save_figure_perhaps(fig, 'ExcWeights-Thresholds_totalKCresp_1odor_contour_reducedSize');


%% Make a 3D plot of iso-surfaces 

[w,a,th] = meshgrid( reducedSet.alphaFactors, reducedSet.weightsFactors, reducedSet.thetaFactors);
reducedSubspace = compareAllParams(1:Excit_cutoff, 1:APL_cutoff, 1:theta_cutoff)/1000;

fig = figure;
ax = axes;
hold on
% for level=(1:2:14)*1e4
for level=(10:10:80)
    isosurface(w,a,th,reducedSubspace, level);
    % isosurface(compareAllParams, level);
end
xlabel('weights increasing by factor')
ylabel('alpha increasing by factor')
zlabel('C_{\theta} increasing')
for ch=[1:length(ax.Children)-2, length(ax.Children)]
    ax.Children(ch).FaceAlpha = 0.4;
end
view(3)

save_figure_perhaps(fig, 'threeParams_totalKCresp_1odor_isosurface_reducedSize');



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

% modify responses of "alphaPrime (alPrime) lobe" by applying excitory-weights and alpha factors
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




%% local functions

function save_figure_perhaps(fig, figureName)
    global saveFigures
    global savepath
    global saveFormats
    if saveFigures
        for format=saveFormats
            saveas(fig,[savepath, figureName, format{:}])
        end
    end
end


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
    ax.Box='off';
    
    
    hold on
    plot(repmat([0.85;1.15],1,numModels), [mean(vanillaResponses,2) , ...
        mean(vanillaDisinh,2)]'/1000 , 'color','#606060')
    plot(repmat([1.85;2.15],1,numModels), [mean(adaptedModelResp,2) , ...
        mean(adaptedModelDisinh,2)]'/1000 , 'color','#606060')
    
    legend({'normal','disinhibited'}, 'fontsize',16)
end

function isequal = isEqualFloat(a,b) 
    % need a dumb little function for float comparison bc matlab fails at 
    % it above
    %isequal = abs(a-b) < 1e-6; 
    isequal = abs(a-b) < eps(a)*5; %even better, x times float accuracy
end