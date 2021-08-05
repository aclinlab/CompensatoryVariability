terpenes = 26:41;
alcohols = 69:86;
esters = 87:110;

numToSample = 15;
numSeqs = 20;
[TerpeneOds, AlcoOds, EstersOds] = deal(cell(1,numSeqs));
for i=1:numSeqs
    TerpeneOds{i} = randsample(terpenes, numToSample, false); % sample without replacement
    AlcoOds{i} = randsample(alcohols, numToSample, false);
    EstersOds{i} = randsample(esters, numToSample, false);
end

save('random15TerpeneAlcoholEsters.mat','TerpeneOds','AlcoOds','EstersOds');