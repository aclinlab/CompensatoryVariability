%% this code imports the already tuned flies, with their input synaptic weights 
%% tuned for activity equalization using: (a) the simple rule (our main model) 
%% (b) the main model (with the extra <xjk>k) (c)equalizing the activities
%% only for the active KCs. This is the data used for panels Fig5F, S3 C&D.
clear all 
clc
CurrentDir=pwd;

FigWin=figure(1);
edges=[0:0.1:4];
% sampling from the experimental data by Turner et al
TurnerExp=  randomTurnerMinis( 12000);
h1=histogram(TurnerExp,edges,'Normalization','probability','DisplayStyle','stairs', 'EdgeColor','k','parent',FigWin);

% the lognormal fit to Turner's data
for i=1:100000
    w(i)=(exp(-0.0507+(0.3527*randn(1))));
end

fig2=figure;
hold on, P_w= histogram(w,edges,'Normalization','probability');
P_lognor=P_w.BinCounts;
close(fig2);

Mu_w=mean(w);
std_w=std(w);

figure(1);
hold on, plot([0:0.1:3.95],(P_lognor./(100000)),'k','LineWidth',1.5);

Ftemp=figure(2);
% import the tuned input weights matrices in the main model
% for each random network, each model fly

for fly=1:20
    
     load(strcat(CurrentDir,'/TunedFlies_allModels',['_fly_wNoise',num2str(fly),num2str(1)]),'thisW_ActivityBasedComp_noxjk')
      W(:,:,fly)=thisW_ActivityBasedComp_noxjk;
      
      % scale the tuned distribution of weights by the ratio of the mean of
      % the lognormal fit for Turner's data: the mean of the tuned distribution
      % this scaling perserves the distribution statistics, it is just used 
      %  so we can compare both distributions in the same domain of W.
      
      w_Eq2= thisW_ActivityBasedComp_noxjk(find(thisW_ActivityBasedComp_noxjk)).* (Mu_w/mean(thisW_ActivityBasedComp_noxjk(find(thisW_ActivityBasedComp_noxjk))) );
      
      h1= histogram(w_Eq2, edges,'LineStyle','-','Normalization','probability');
      countsEq2(:,fly) =h1.BinCounts./size(w_Eq2,1);
   
end

countsEq2_avg_flies=mean(countsEq2,2);
close(Ftemp);
figure(1);
hold on, plot([0:0.1:3.95],countsEq2_avg_flies,'color','b','LineWidth',1.5);

%% for FigS3 C&D
% open a new window to plot the Turner' data and the lognormal fit. 

FigWin2=figure(3);
h1=histogram(TurnerExp,edges,'Normalization','probability','DisplayStyle','stairs', 'EdgeColor','k','parent',FigWin2);
hold on, plot([0:0.1:3.95],(P_lognor./(100000)),'k','LineWidth',1.5);
hold off

FTemp=figure(2)
% do the same for the alternate rules for tuning W: FigS3 C&D 
for fly=1:20
    
     load(strcat(CurrentDir,'/TunedFlies_allModels',['_fly_wNoise',num2str(fly),num2str(1)]),'thisW_ActivityBasedComp','thisW_ActivityBasedComp_wHy')     
     W_1(:,:,fly)=thisW_ActivityBasedComp; 
     W_2(:,:,fly)=thisW_ActivityBasedComp_wHy;
     
     % scale the tuned distribution of weights by the ratio of the mean of
     % the lognormal fit for Turner's data: the mean of the tuned distribution
     w_Eq2_1= thisW_ActivityBasedComp(find(thisW_ActivityBasedComp)).* (Mu_w/mean(thisW_ActivityBasedComp(find(thisW_ActivityBasedComp))) );
     w_Eq2_2= thisW_ActivityBasedComp_wHy(find(thisW_ActivityBasedComp_wHy)).* (Mu_w/mean(thisW_ActivityBasedComp_wHy(find(thisW_ActivityBasedComp_wHy))) );
     
     h11= histogram(w_Eq2_1, edges,'LineStyle','-','Normalization','probability');
     countsEq2_1(:,fly) =h11.BinCounts./size(w_Eq2_1,1);
     
     h12= histogram(w_Eq2_2, edges,'LineStyle','-','Normalization','probability');
     countsEq2_2(:,fly) =h12.BinCounts./size(w_Eq2_2,1);
   
end

countsEq21_avg_flies=mean(countsEq2_1,2);
countsEq22_avg_flies=mean(countsEq2_2,2);

close(FTemp);

figure(3);
hold on, plot([0:0.1:3.95],countsEq21_avg_flies,'color',[0/255 0/255  139/255],'LineWidth',1.5);
hold on, plot([0:0.1:3.95],countsEq22_avg_flies,'color',[102/255 178/255  255/255],'LineWidth',1.5);

