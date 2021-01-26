%% this code imports the stored tuned flies from the current directory, 
%% it then plots the distribution of the KCs spiking thresholds
%% tuned for activity equalization or equalizing firing probabilities. 
%% this is the data used in panels Fig5F and FigS4C

clear all 
clc
CurrentDir=pwd;

FigWin=figure(1);

% spiking thresholds limits in Turner et al. data
edges=[0:2.5:40];

TurnerExp_theta= randomTurnerSpikingThresh( 2000);

histogram(TurnerExp_theta,edges,'Normalization','probability','DisplayStyle','stairs','Parent',FigWin,'LineWidth',1.5,'EdgeColor','k');

% the Gaussian fit to the experimental data
for i=1:2000
    theta(i)=(21.5+5.6.*randn(1));
end

fig2=figure;
P_theta= histogram(theta,edges,'Normalization','probability');
P_theta_G= P_theta.BinCounts;
close(fig2);

Mu_T=mean(theta);
std_T=std(theta);

figure(1);
hold on, plot([0:2.5:37.5],(P_theta_G./(2000)),'k','LineWidth',1.5);


fig3=figure;
% import the tuned spiking thresholds in our main model and from the
% Kennedy's inspired model(firing probabilities equalization)
% for each random network, each model fly

for fly=1:20
        load( strcat(CurrentDir,'/TunedFlies_allModels',['_fly_wNoise',num2str(fly),num2str(1)]));
        Theta(:,fly)=theta_Activity_homeo; 
        T= theta_Activity_homeo.* (Mu_T/mean(theta_Activity_homeo) );
        h1= histogram(T, edges,'LineStyle','-','Normalization','probability');
        counts_theta(:,fly) =h1.BinCounts./(2000);
        
        T_Kenn= theta_Kenn.* (Mu_T/mean(theta_Kenn) );
        h1_Kenn= histogram(T_Kenn, edges,'LineStyle','-','Normalization','probability');
        countsKenn(:,fly) =h1_Kenn.BinCounts./(2000);
        
        

end

counts_theta_avg_flies=mean(counts_theta,2);
countsKenn_avg_flies=mean(countsKenn,2);
 

figure(1);
plot([0:2.5:37.5],counts_theta_avg_flies,'color','m','LineWidth',1.5);
plot([0:2.5:37.5],countsKenn_avg_flies,'color',[139/255 0/255 139/255],'LineWidth',1.5);


