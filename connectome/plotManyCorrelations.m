allParams = [numPNsPerKC, numPNKCsynapsesPerKC, ...
    meanWPerKC, PNDistToPedPerKC, numAPLKCsynapsesPerKC];

numParams = size(allParams,2);
numKCtypes = length(indicesByTypeWithinAllKCindices);

% cell arrays to store the results for each KC subtype
correlations = cell(numKCtypes,1);
pvalues = zeros(numParams,numParams,numKCtypes);

% loop through KC subtypes
for i=1:numKCtypes
    [correlations{i}, pvalues(:,:,i)] = corr(allParams(indicesByTypeWithinAllKCindices{i},:));
end

% apply Holm-Bonferroni correction to p-values
% separate out upper triangle and lower triangle because they are not
% actually separate comparisons (i.e. for n x n matrix of pairwise
% comparisons, there are only n*(n-1)/2 multiple comparisons)
[pvaluestriu, pvaluestril] = deal(zeros(size(pvalues)));
for i=1:numKCtypes
    pvaluestriu(:,:,i) = triu(pvalues(:,:,i),1);
    pvaluestril(:,:,i) = tril(pvalues(:,:,i),-1);
end

[pvaluesHBtriu, ~] = bonf_holm(pvaluestriu);
[pvaluesHBtril, ~] = bonf_holm(pvaluestril);

pvaluesHB = pvaluesHBtriu + pvaluesHBtril;
% now need to NaN the diagonal
for i=1:numKCtypes
    x = pvaluesHB(:,:,i);
    x(eye(size(x))==1) = NaN;
    pvaluesHB(:,:,i) = x;
end

figure
% loop through KC subtypes
for i=1:numKCtypes
    subplot(1,numKCtypes,i)
    imagesc(correlations{i}, [-1 1]);
    title(types{i});
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',[],'YTick',[])
    axis equal, axis tight
    hold on
    for j=1:size(pvalues(:,:,i),1)
        for k=1:size(pvalues(:,:,i),2)
            if (pvaluesHB(j,k,i) < 0.05) 
                % p-value with Holm-Bonferroni correction for multiple comparisons
                % put a black dot if the correlation is significant
                scatter(j,k,'k','filled');
            end
        end
    end
end
colormap(redblue)

function result = upperRight(input)
if size(input,1)~=size(input,2)
    error('first 2 dimensions need to be square (equal dims)');
end
if ndims(input)==3
    result = zeros(size(input,1)^2-size(input,1),size(input,3));
    for i=1:size(input,3)
        thisSlice = input(:,:,i);
        result(:,i) = thisSlice(~eye(size(input,1)));
    end
else
    result = input(~eye(size(input)));
end
end