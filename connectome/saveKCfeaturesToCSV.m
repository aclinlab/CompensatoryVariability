paramNames = {'numPNsPerKC', 'numPNKCsynapsesPerKC', ...
    'meanWPerKC', 'PNDistToPedPerKC', 'numAPLKCsynapsesPerKC'};
table = array2table(allParams,'VariableNames', paramNames);
table.KCindex = allKCindices;
table.bodyId = KCskel.getBodyIDs(allKCindices);
table.subtype = KCskel.getTypes(allKCindices);
writetable(table,'Fig7_connectome_data.csv');