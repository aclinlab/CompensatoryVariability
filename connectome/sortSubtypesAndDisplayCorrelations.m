%% pull out KC indices by subtype

allowedTypes = {
    'KCg-m_R';
    'KCab-s_R';
    'KCab-m_R';
    'KCab-c_R';
    'KCa''b''-ap2_R';
    'KCa''b''-m_R';
    };

[indices, bodyids] = KCskel.getKCsOfTypes(allowedTypes);
%indices is a cell array
types = allowedTypes;
indicesByType = indices;
% double check we throw away the excluded KCs
KCsToExclude = csvread('KCsToExclude.csv');
for i=1:length(indicesByType)
    indicesByType{i} = setdiff(indicesByType{i},KCsToExclude);
end
indicesByTypeWithinAllKCindices = cell(length(indicesByType),1);
for i=1:length(indicesByType)
    % [C,IA,IB] = intersect(A,B) returns index vectors IA and IB such
    % that C = A(IA) and C = B(IB).
    [~,indicesByTypeWithinAllKCindices{i},~] = intersect(allKCindices, indicesByType{i});
end

% %% display correlations
% 
% 
disp('w vs N');
plotSubtypeCorrelations(numPNsPerKC, '# PNs per KC', ...
    meanWPerKC, '# synapses per PN-KC connection', ...
    indicesByTypeWithinAllKCindices, types, 'violins')
saveas(gcf,'wvsN.svg','svg')
% 
disp('dist to post-ped vs N')
plotSubtypeCorrelations(numPNsPerKC, '# PNs per KC', ...
    PNDistToPedPerKC, 'dist. to post. ped. (um)', ...
    indicesByTypeWithinAllKCindices, types, 'violins')
saveas(gcf,'dist-to-post-ped-vs-N.svg','svg')

disp('alpha vs all PN-KC synapses')
plotSubtypeCorrelations(numPNKCsynapsesPerKC, '# PN-KC synapses per KC', ...
    numAPLKCsynapsesPerKC, '# APL-KC synapses per KC', ...
    indicesByTypeWithinAllKCindices, types, 'scatter')
saveas(gcf,'alpha-vs-totalPNKCsyns.svg','svg')

disp('all PN-KC synapses vs N')
plotSubtypeCorrelations(numPNsPerKC, '# PNs per KC',...
    numPNKCsynapsesPerKC, '# PN-KC synapses per KC', ...
    indicesByTypeWithinAllKCindices, types, 'violins')
saveas(gcf,'totalPNKCsyns-vs-N.svg','svg')

