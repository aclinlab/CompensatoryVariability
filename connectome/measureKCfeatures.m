% %%
allowedTypes = {
    'KCg-m_R';
    'KCab-s_R';
    'KCab-m_R';
    'KCab-c_R';
    'KCa''b''-ap2_R';
    'KCa''b''-m_R';
    };

[indices, bodyids] = KCskel.getKCsOfTypes(allowedTypes);
KCsToExclude = csvread('KCsToExclude.csv');
allKCindices = sort(setdiff(cell2mat(indices),KCsToExclude));
numKCs = length(allKCindices);

%% allocate memory
% these variables will store the measurements for all KCs. 
%
numAPLKCsynapsesPerKC = zeros(numKCs,1); % number of APL->KC synapses per KC
numPNsPerKC = zeros(numKCs,1); % number of projection neurons that connect to each KC (ignoring connections of <=2 synapses)
numPNKCsynapsesPerKC = zeros(numKCs,1); % total number of PN-KC synapses per KC
numKCAPLsynapsesPerKC = zeros(numKCs,1); % number of APL->KC synapses per KC

% distance of each PN synapse to the posterior boundary of the peduncle, 
% first averaged across all synapses for each PN, then averaged across all
% PNs for that KC
PNDistToPedPerKC = zeros(numKCs,1); 
% as above but not averaged across all PNs for each KC - preserve
% individual PNs 
PNDistToPedPerKCPerPN = cell(numKCs,1);

% individual weights per PN-KC connection
wPerPNKCconnection = cell(numKCs,1);

%% load KCs

f = waitbar(0,'Loading KC .mat files');

for i=1:numKCs
    if mod(i,100)==0
        waitbar(i/numKCs, f);
    end
    load(strcat('kc',num2str(allKCindices(i)),'.mat'));
    kcs(i) = obj;
end
    
close(f)

%% extract data about each KC
    
% loop through all selected KCs
f = waitbar(0,'Measuring KC features');
for i=1:numKCs
    if mod(i,20)==0
        waitbar(i/numKCs, f);
    end
    
    % get # APL-KC synapses in the calyx
    
    APLtoKC = find(kcs(i).synSet(2).uniquePartners == 425790257); %425790257 is the bodyId of APL
    % use input synapses that are:
    % 1. from APL
    % 2. posterior to the 160 um mark
    % 3. not in the peduncle (i.e. the PED(R) ROI)
    % we use this strange definition of "calyx" because some calyx APL-KC
    % synapses are not annotated as being in the CA(R) ROI
    APLtoKCsyns = find((kcs(i).synSet(2).uniquePartnerIndices == APLtoKC) &...
        (kcs(i).synSet(2).synLocs(kcs(i).synSet(2).synLocsOfOrigSyn,2)<(160/.008)) &...
        ~contains(kcs(i).synSet(2).rawData.upstream_roi,'PED(R)'));
    numAPLKCsynapsesPerKC(i) = length(APLtoKCsyns);

    KCtoAPL = find(kcs(i).synSet(1).uniquePartners == 425790257); %425790257 is the bodyId of APL
    % use input synapses that are:
    % 1. from APL
    % 2. posterior to the 160 um mark
    % 3. not in the peduncle (i.e. the PED(R) ROI)
    % we use this strange definition of "calyx" because some calyx APL-KC
    % synapses are not annotated as being in the CA(R) ROI
    KCtoAPLsyns = find((kcs(i).synSet(1).uniquePartnerIndices == KCtoAPL) &...
        (kcs(i).synSet(1).synLocs(kcs(i).synSet(1).synLocsOfOrigSyn,2)<(160/.008)) &...
        ~contains(kcs(i).synSet(1).rawData.upstream_roi,'PED(R)'));

    % not corrected for confidence
    numKCAPLsynapsesPerKC(i) = length(KCtoAPLsyns);
    
    % count PN inputs
    inputPNsForThisKC = kcs(i).inputPNs();
    numPNsPerKC(i) = length(inputPNsForThisKC);

    % count total PN-KC synapses
    synapsesPerPN = zeros(numPNsPerKC(i),1);
    for j=1:length(synapsesPerPN)
        PNtoKCsyns = find(kcs(i).synSet(2).uniquePartnerIndices == inputPNsForThisKC(j));
        
        synapsesPerPN(j) = length(PNtoKCsyns);
        
    end
    numPNKCsynapsesPerKC(i) = sum(synapsesPerPN);

    % measure distance from PNs to posterior boundary of peduncle
    distancesToPostPed = kcs(i).meanPNdistByPN('real','postPed')*0.008; % 0.008 um per pixel
    PNDistToPedPerKCPerPN{i} = distancesToPostPed(isfinite(distancesToPostPed));
    PNDistToPedPerKC(i) = mean(distancesToPostPed(isfinite(distancesToPostPed)));
        
    wPerPNforThisKC = kcs(i).numSynapsesPerPN;
    wPerPNKCconnection{i} = wPerPNforThisKC(isfinite(distancesToPostPed)); % ensure same # of PNs to allow correlations
    
end
close(f); % waitbar
meanWPerKC = numPNKCsynapsesPerKC./numPNsPerKC;


