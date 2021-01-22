function slopes = plotSubtypeCorrelations(x, xtitle, y, ytitle, indices, types, graphtype)
% plotSubtypeCorrelations   
%   plot correlations of 2 KC parameters for each subtype of KC
%   plotSubtypeCorrelations(x, xtitle, y, ytitle, indices, types, graphtype)
%   x = array (numKCs x 1) containing one parameter measured for each KC
%   xtitle = a string used to label the x-axis in the resulting plots,
%       describing the parameter in x
%   y = array (numKCs x 1) containing another parameter measured for each KC
%   ytitle = a string used to label the y-axis in the resulting plots,
%       describing the parameter in y
%   indices = a cell array (3 x 1) where each cell contains an array with
%       the indices within x and y of each subtype of KC. e.g. indices{1}
%       might be [1 3 4 5 7 8 10 11 12...]. If calling this function after
%       running measureKCfeatures, this is probably indicesByTypeWithinAllKCindices
%   types = a cell array of strings labeling the subtypes of KCs indexed in
%       indices. If calling this function after running measureKCfeatures,
%       this is the variable types
%   graphtype = a string with only 2 valid options:
%       'violins' - graph violin plots. Use this where one parameter has
%          only a few possible integer values, e.g. numPNsPerKC
%       'scatter' - scatter plot. Use this where both parameters can have a
%          a wide range of values
%   Example call: plotSubtypeCorrelations(numPNsPerKC, '# PNs per KC', ...
%     meanWPerKC, '# synapses per PN-KC connection', ...
%     indicesByTypeWithinAllKCindices, types, 'violins')


color = 'r';
figure; % draw empty figure
slopes = zeros(length(indices),1);
for i=1:length(indices)
     subplot(1,length(indices),i) % the graphs will appear side by side
% .    figure
    % draw the graph: violin plot or scatter plot
    switch graphtype
        case 'violins'
            distributionPlot(y(indices{i}),...
                'groups',x(indices{i}),'globalNorm',1,'showMM',3)
            set(gca,'XTick',1:max(x(indices{i})));
        case 'scatter'
            scatter(x(indices{i}), y(indices{i}), 6, 'k','MarkerEdgeAlpha', 0.2);
        otherwise
            error('unknown type');
    end
    
    % get a linear fit
    p = polyfit(x(indices{i}), y(indices{i}) ,1);
    p_fix = polyfix(x(indices{i}), y(indices{i}) ,1,0,0); % fix through the origin
    hold on
    slopes(i) = p(1)/mean(y(indices{i}));
    
    % draw the linear fit
    x_line = min(x(indices{i})):max(x(indices{i}));
    plot(x_line,p(1)*x_line+p(2),'Color',color,'LineWidth',1)
    % uncomment this line to draw the trend line going through the origin
    % plot(x_line,p_fix(1)*x_line+p_fix(2),'Color','k','LineWidth',0.5)
    
    
    % format the graph
    ylim([0 max(y(indices{i}))]);
    set(gca,'TickLength',[.03 .03],'TickDir','out')
    set(gca,'FontSize',6)
    set(gca,'units','centimeters','position',[1+2.5*(i-1),1,2.1,2.1])
    xlabel(xtitle,'FontSize',6);
    if i==1
        ylabel(ytitle,'FontSize',6);
    end
    title(types{i});
    
    % print out the correlation coefficient and the p-value
    [r,pvalue] = corrcoef(x(indices{i}), y(indices{i}));
    fprintf('Pearson correlation value for %s: %.2f\n', types{i}, r(2,1));
    fprintf('p-value of correlation for %s: %d\n', types{i}, pvalue(2,1)); 
    % note this is not corrected for multiple comparisons - see
    % plotManyCorrelations for Holm-Bonferroni correction
    
end

end
