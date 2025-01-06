function [PNactivity, spontanFRate, odorList] = propagateORNs2PNs(varargin)
%% get PN activity from ORNs. Accepts ORN input as argument, as a K-by-M 
% format with K response patterns in M ORNs (or glomeruli)
% When no argument is given, it will operate on the Hallem dataset

if length(varargin)==0
    tb = readtable('hallem and carlson 2006.xlsx');
    tb = renamevars(tb, 'Var1','odorName');
    tb = removevars(tb, 'Var2');
    
    tb(cellfun(@isempty,tb.odorName),:) = []; %remove empty line(s) (expect one)
    spontanFRate = tb(end,:);
    
    tb(end,:)=[];
    odorList = tb.odorName;
    ORNresponses = table2array(removevars(tb, 'odorName'));
    spontanFRate = table2array(spontanFRate(1,2:end) );
    % floor the negative deltas at 0 (instead of adding to the spont. rate)
    % complex debate around this
    ORNresponses(ORNresponses<0) = 0; 
    % or instead add the spontaneous firing rate to the deltas
    % ORNresponses = ORNresponses + repmat(spontanFRate, size(ORNresponses,1),1);
    % a third option is implemented in the calculation, and is probably the
    % correct one
else
    ORNresponses = varargin;
end


m = 10.63; %gain of lateral inhibition in AL
Rmax = 165; %maximum PN response
sigma = 12; %non-linearity parameter of ORN to PN response function
s = m * sum(ORNresponses,2)/190;

% third option of how to deal with deltas data: set PN response to 0 if the
% presynaptic ORN is suppressed from baseline, but keep suppressed ORNs for
% calculating the lateral inhibition. This makes sense especially 
% PNactivity = double(ORNresponses>0) .* Rmax .* ORNresponses.^1.5 ./(ORNresponses.^1.5 +s.^1.5 +sigma.^1.5);
PNactivity = Rmax .* ORNresponses.^1.5 ./(ORNresponses.^1.5 +s.^1.5 +sigma.^1.5);


end