classdef KCskel < neurSkel
    properties
        postPedNeighbors
        realDistBetweenPNInputs
        EDBetweenPNInputs
        realDistPostPedToPNInputs
        EDPostPedToPNInputs
    end
    
    properties (Constant)
        outputIndex = 1
        inputIndex = 2
        PNIndex = 3
        postPedIndex = 5 % for legacy reasons
    end
    
    % Follow this convention for synSet:
    % 1: KC output synapses, i.e. KC->X
    % 2: KC input synapses, i.e. X->KC
    % 3: PN-KC synapses, i.e. PN->KC (subset of 2, create using extractPNsynapses
    % 5: calyx-peduncle boundary according to hemibrain ROIs
    
    % dimension conventions:
    % in the connectome dimensions are like 1,2,3
    % 1 is lateral->medial
    % 2 is posterior->anterior
    % 3 is dorsal->ventral
    % in Matlab when plotting you have to swap x and y
    % so in drawSkeleton they are plotted like 2,1,3 
    
    methods
        
        function obj = KCskel(filename)
            obj = obj@neurSkel(string(filename));
        end
        
        function [] = saveKCskel(obj, varargin)
            if nargin>1
                % user gives file name to save the object as
                filename = varargin{1};
            else
                % extract the number of the KCskel.csv file
                % use that to automatically generate a name for the .mat
                % file
                filename = strcat('kc',num2str(obj.getKCindex),'.mat');
            end
            save(filename, 'obj');
            
        end
        
        function number = getKCindex(obj)
            % extract the number of the KCskel.csv file stored in the
            % object
            [startIndex, endIndex] = regexp(obj.filepath,'KCskel[0-9]{1,4}\.csv', 'once'); % see https://en.wikipedia.org/wiki/Regular_expression
            if ~isempty(startIndex)
                charfilepath = char(obj.filepath);
                number = str2double(charfilepath((startIndex+6):(endIndex-4)));
            else
                error('Could not extract KC number');
            end
        end
        
        function result = inputPNs(obj)
            % return a list of input PNs (INDEXED WITHIN
            % obj.synSet(2).uniquePartners) that form more than 2 synapses
            % on this KC
            PNs = readtable('monoglom_PNs_v1.1.csv');
            [~,result,~] = intersect(obj.synSet(2).uniquePartners, PNs{:,'bodyId'});
            for i=1:length(result)
                if (nnz(obj.synSet(2).uniquePartnerIndices == result(i)) <= 2)
                    result(i) = nan;
                end
            end
            result(isnan(result)) = [];
        end
        
        function obj = definePosteriorPeduncleEdge(obj)
            % creates a new synSet structure at index synSetIndex
            % (convention is this should be 5)
            % at the more anterior node of the link that straddles the
            % posterior boundary of the peduncle
            % relies on some assumptions, might run into some edge
            % cases....
            peduncleIndex = find(strcmp(obj.ROIskey,'PED(R)'));
            linksInROI = obj.linkROIs==peduncleIndex;
            nodesInROI = obj.nodeROIs==peduncleIndex;
            if ~isempty(linksInROI)
                % use obj. links to index into nodesInROI, which is 0 for
                % outside ROI, 1 for inside ROI. This produces a numNodes
                % x 2 logical matrix
                % ignore links that have 0 values
                linkIndices = 1:size(obj.links,1);
                nonzeroLinkIndices = linkIndices(all(obj.links,2));
                
                endpointsInROI = nodesInROI(obj.links(nonzeroLinkIndices,:));
                % find the links where one endpoint is in the ROI and the
                % other one is not
                inOutLinks = nonzeroLinkIndices((endpointsInROI(:,1) - endpointsInROI(:,2))~=0);
                if length(inOutLinks)~=2
                    disp(strcat('warning: length(inOutLinks)=', num2str(length(inOutLinks)),', ', obj.filepath, ' might not be entirely in the peduncle ROI'))
                end
                
                linkMidpoints = zeros(length(inOutLinks), 3);
                % get the most posterior link
                for i=1:length(inOutLinks)
                    linkMidpoints(i,:) = mean([obj.nodes(obj.links(inOutLinks(i),1)).coords; ...
                        obj.nodes(obj.links(inOutLinks(i),2)).coords], 1);
                end
                antPostDim = linkMidpoints(:,2);
                [~,mostPostInOutLink] = min(antPostDim);
                calyxPeduncleLink = inOutLinks(mostPostInOutLink);
                % get the more anterior node from the most posterior link
                % (this is the most posterior node of the peduncle)
                calyxPeduncleNodes = cell2mat({obj.nodes(obj.links(calyxPeduncleLink,:)).coords}');
                [~,mostAntNodeIndex] = min(calyxPeduncleNodes(:,2));
                mostAntNode = obj.links(calyxPeduncleLink,mostAntNodeIndex);
                obj.synSet(obj.postPedIndex).synLocs = obj.nodes(mostAntNode).coords;
                obj.synSet(obj.postPedIndex).synLinks = calyxPeduncleLink;
                [~,rootname,~] = fileparts(obj.filepath);
                obj.synSetNames{obj.postPedIndex} = strcat(rootname,'-posteriorPeduncleEdge');
            else
                error('No links in the ROI not defined')
            end
        end
        
        function obj = extractPNsynapses(obj)
            % create a synSet just for PN inputs (goes into array of synSet
            % structures, called obj.synSet - into index synSetIndex)
            
            % keep only the uniquePartners that are PNs
            inputPNindices = obj.inputPNs();
            obj.synSet(obj.PNIndex).uniquePartners = obj.synSet(2).uniquePartners(inputPNindices);
            
            % keep only the synapses that come from PNs
            PNpartnerIndices = ismember(obj.synSet(2).uniquePartnerIndices, inputPNindices);
            obj.synSet(obj.PNIndex).uniquePartnerIndices = obj.synSet(2).uniquePartnerIndices(PNpartnerIndices);
            % change uniquePartnerIndices to be e.g. 1 2 3 4 not 56 89 100 153 etc
            [~,~,obj.synSet(obj.PNIndex).uniquePartnerIndices] = unique(obj.synSet(obj.PNIndex).uniquePartnerIndices);
            
            obj.synSet(obj.PNIndex).rawData = obj.synSet(2).rawData(PNpartnerIndices,:); % Note this doesn't change the indexing in column 'Var1'

            obj.synSet(obj.PNIndex).synLocsOfOrigSyn = obj.synSet(2).synLocsOfOrigSyn(PNpartnerIndices);
            [PNsynLocs,~,ic] = unique(obj.synSet(obj.PNIndex).synLocsOfOrigSyn);
            % change synLocsOfOrigSyn to be e.g. 1 2 3 4 not 76 78 79 etc
            %prevSynIndex = obj.synSet(synSetIndex).synLocsOfOrigSyn;
            obj.synSet(obj.PNIndex).synLocsOfOrigSyn = ic;
            %PNsynLocs = unique(obj.synSet(synSetIndex).synLocsOfOrigSyn);
            
            % synGroupings indexing is incorrect for the moment. Fix it
            % later if need be
%             obj.synSet(synSetIndex).synGroupings = obj.synSet(2).synGroupings(PNsynLocs);
%             for i=1:length(obj.synSet(synSetIndex).synGroupings)
%                 obj.synSet(synSetIndex).synGroupings{i} = PNsynLocs(obj.synSet(synSetIndex).synGroupings{i});
%             end
            obj.synSet(obj.PNIndex).synLocs = obj.synSet(2).synLocs(PNsynLocs,:);
            obj.synSet(obj.PNIndex).synLinks = obj.synSet(2).synLinks(PNsynLocs);
            [~,rootname,~] = fileparts(obj.filepath);
            obj.synSetNames{obj.PNIndex} = strcat(rootname,'-PNinputs');
        end
            
        function obj = getDistBetweenPNInputs(obj)
            % isolates all PN input synapses and saves the real and
            % electrotonic distances between them
            
            % get PNs
            obj = obj.extractPNsynapses;
            obj = obj.getSynSetGraph(obj.PNIndex);
            
            % get all distances between PN inputs
            obj.EDBetweenPNInputs = obj.getAllSynDistWithinSynSet(obj.PNIndex, 1e7, 3);
            obj.realDistBetweenPNInputs = obj.getAllSynDistWithinSynSet(obj.PNIndex, 1e7, 2);
        end
        
        function obj = getDistPostPedToPNInputs(obj)
            % get the nearest neighbor PN input synapses to the posterior
            % boundary of the peduncle
            obj.postPedNeighbors = obj.getNeighborsOnOtherGraph([obj.postPedIndex obj.PNIndex]);
            
            % get distances from the peduncle to PN inputs
            obj.EDPostPedToPNInputs = obj.getDistBetweenSynSets(obj.EDBetweenPNInputs, obj.postPedNeighbors, 1);
            obj.realDistPostPedToPNInputs = obj.getDistBetweenSynSets(obj.realDistBetweenPNInputs, obj.postPedNeighbors, 0);
        end

        function result = meanPNdistByPN(obj, varargin)
            % optional argument is whether to use electrotonic or real distance
            % example call:
            % kc.meanPNdistByPN('real')
            % default: electrotonic distance
            if nargin>1
                switch varargin{1}
                    case 'electrotonic'
                        distances = obj.EDPostPedToPNInputs;
                    case 'real'
                        distances = obj.realDistPostPedToPNInputs;
                    otherwise
                        error('unrecognised input to function');
                end
            else
                distances = obj.EDPostPedToPNInputs;
            end
            
            numPNs = length(obj.synSet(3).uniquePartners);
            result = zeros(numPNs,1);
            for i=1:numPNs
                % find the synapses that have this PN as the partner
                % then use that logical indexing to find the matching
                % synLocs in synLocsOfOrigSyn
                % then use synLoc indexes to index in distances
                thisPNdistances = distances(obj.synSet(3).synLocsOfOrigSyn(obj.synSet(3).uniquePartnerIndices==i));
                result(i) = mean(thisPNdistances(isfinite(thisPNdistances))); % use isfinite to ignore NaN and Inf
            end
        end
        
        function result = numSynapsesPerPN(obj)
          
            numPNs = length(obj.synSet(3).uniquePartners);
            result = zeros(numPNs,1);
            for i=1:numPNs
                result(i) = nnz(obj.synSet(3).uniquePartnerIndices==i);
            end
                        
        end
                
        function [axHandle] = drawSkeleton(obj, varargin)
            % This is for real neurons but beware of very large skeletons
            
            if nargin>1
                ax = varargin{1};
                if nargin>2
                    neuriteColor = varargin{2};
                else
                    neuriteColor = 'k';
                end
            else
                figure;
                ax = gca;
                neuriteColor = 'k';
            end
            
            % ignore links that have 0 values
            linkIndices = 1:size(obj.links,1);
            nonzeroLinkIndices = linkIndices(all(obj.links,2));
            if length(nonzeroLinkIndices)>10000
                resp = input('This is a huge skeleton (>10000 nodes). Are you sure you want to continue? press n to exit');
                if (resp=='n')
                    return;
                end
            end
            
%             title(obj.filepath);
            %figHandle = figure();
            hold on;
            % draw skeleton links
            for j=nonzeroLinkIndices
                coords(1,:) = obj.nodes(obj.links(j,1)).coords; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).coords; % 1x3
                % coords is now a 2x3 array
                peduncleIndex = find(strcmp(obj.ROIskey,'PED(R)'));
%                 if ~isempty(obj.linkROIs) && obj.linkROIs(j)==peduncleIndex
%                     % draw the peduncle in blue
%                     color = 'b';
%                 else
                    color = neuriteColor;% other segments are black
%                 end
                line(ax,coords(:,2)',coords(:,1)',-coords(:,3)','Color',color,'LineWidth',mean(obj.skel.radius(obj.links(j,:)))/20);
            end
            cmap = colormap;
            %colorList = ['m','c','g','r','b','k','y'];
            inputPNs = obj.inputPNs();
            

            % draw input PNs
            for i=1:length(inputPNs)
                synapsesForThisPN = find(obj.synSet(2).uniquePartnerIndices == inputPNs(i));
                for j=1:length(synapsesForThisPN)
                    color = cmap(ceil(i/length(inputPNs)*size(cmap,1)),:);
                    % this commented-out code is to draw the locations of
                    % the PN pre-synaptic terminals
%                     x = obj.synSet(2).rawData{obj.synSet(2).synGroupings{synapsesForThisPN(j)}(1),'uy'}; %(synapsesForThisPN(j),2);
%                     y = obj.synSet(2).rawData{obj.synSet(2).synGroupings{synapsesForThisPN(j)}(1),'ux'};
%                     z = -obj.synSet(2).rawData{obj.synSet(2).synGroupings{synapsesForThisPN(j)}(1),'uz'};
                    x = obj.synSet(2).synLocs(obj.synSet(2).synLocsOfOrigSyn(synapsesForThisPN(j)),2);
                    y = obj.synSet(2).synLocs(obj.synSet(2).synLocsOfOrigSyn(synapsesForThisPN(j)),1);
                    z = -obj.synSet(2).synLocs(obj.synSet(2).synLocsOfOrigSyn(synapsesForThisPN(j)),3);
                    plot3(ax,x,y,z,'o','MarkerSize',4,'MarkerEdgeColor',color,'MarkerFaceColor',color); % need to change this when having multiple synSets!
                    
                end
            end
            
            % draw APL synapses
            APLtoKC = find(obj.synSet(2).uniquePartners == 425790257);
            APLtoKCsyns = find(obj.synSet(2).uniquePartnerIndices == APLtoKC);
            for j=1:length(APLtoKCsyns)
                x = obj.synSet(2).synLocs(obj.synSet(2).synLocsOfOrigSyn(APLtoKCsyns(j)),2);
                y = obj.synSet(2).synLocs(obj.synSet(2).synLocsOfOrigSyn(APLtoKCsyns(j)),1);
                z = -obj.synSet(2).synLocs(obj.synSet(2).synLocsOfOrigSyn(APLtoKCsyns(j)),3);
                plot3(ax,x,y,z,'o','MarkerSize',6,'MarkerEdgeColor','r');
            end
%             KCtoAPL = find(obj.synSet(1).uniquePartners == 425790257);
%             KCtoAPLsyns = find(obj.synSet(1).uniquePartnerIndices == KCtoAPL);
%             for j=1:length(KCtoAPLsyns)
%                 x = obj.synSet(1).synLocs(obj.synSet(1).synLocsOfOrigSyn(KCtoAPLsyns(j)),2);
%                 y = obj.synSet(1).synLocs(obj.synSet(1).synLocsOfOrigSyn(KCtoAPLsyns(j)),1);
%                 z = -obj.synSet(1).synLocs(obj.synSet(1).synLocsOfOrigSyn(KCtoAPLsyns(j)),3);
%                 plot3(ax,x,y,z,'x','MarkerSize',6,'MarkerEdgeColor','g');
%             end
            
            % all synapses - commented out
%             for i=1:2
%                 %color = colorList(i); %cmap(i*32,:);
%                 for j=1:length(obj.synSet(i).synLinks)
%                     x = obj.synSet(i).synLocs(j,2);
%                     y = obj.synSet(i).synLocs(j,1);
%                     z = -obj.synSet(i).synLocs(j,3);
%                     plot3(ax,x,y,z,'o','MarkerSize',2,'MarkerEdgeColor','m'); % need to change this when having multiple synSets!
%                 end
%             end
%             
% %             % output synapses that are not KCtoAPL synapses
% %             % ugh... can't figure this out... come back to it.
% %             for j=1:length(obj.synSet(1).synLinks)
% %                 x = obj.synSet(1).synLocs(j,2);
% %                 y = obj.synSet(1).synLocs(j,1);
% %                 z = -obj.synSet(1).synLocs(j,3);
% %                 plot3(x,y,z,'o','MarkerSize',2,'MarkerEdgeColor','m'); % need to change this when having multiple synSets!
% %             end
            
            % 
            
            % draw posterior edge of peduncle
            if length(obj.synSet)>=5
                if size(obj.synSet(5).synLinks==1)
                    x = obj.synSet(5).synLocs(2);
                    y = obj.synSet(5).synLocs(1);
                    z = -obj.synSet(5).synLocs(3);
                    plot3(ax,x,y,z,'square','MarkerSize',6,'MarkerEdgeColor','m');
                end
            end
%             
            
%             linksPerNode = obj.linksPerNode();
%             nodesWithBranches = find(linksPerNode>2);
%             
%             for i=1:length(nodesWithBranches)
%                 loc = obj.nodes(nodesWithBranches(i)).coords;
%                 plot3(loc(2),loc(1),-loc(3),'o','Color','k');
%             end
            axis equal, axis tight
            axHandle = ax;
            view(24,33);
        end

    end
    
    methods (Static)
        function [] = calyxAxisLimits()
            ylim([9875 20000])
            xlim([9875 16268])
            zlim([-16100 -7900])
        end
        
        function ids = getBodyIDs(indices)
            KCs = csvread('KCids_v1.1.csv');
            ids = KCs(indices+1,1); %off-by-one because Python indexes from 0 and Matlab from 1
        end
        
        function indices = getIndices(bodyIds)
            KCs = csvread('KCids_v1.1.csv');
            [~,ia,~] = intersect(KCs,bodyIds);
            indices = ia-1; %off-by-one because Python indexes from 0 and Matlab from 1
        end
        
        function types = getTypes(indices)
            KCids = KCskel.getBodyIDs(indices);
            types = cell(length(KCids),1);
            allNeurons = readtable('traced-neurons-v1.1.csv');
            for i=1:length(KCids)
                types(i) = allNeurons.instance(allNeurons.bodyId==KCids(i));
            end
        end
        
        function [indices, bodyIds] = getKCsOfTypes(types)
            % types must be a cell array of strings
            allowedTypes = {
                'KCg-m_R';
                'KCg-t_R';
                'KCg-d_R';
                'KCab-p_R';
                'KCab-s_R';
                'KCab-m_R';
                'KCab-c_R';
                'KCa''b''-ap1_R';
                'KCa''b''-ap2_R';
                'KCa''b''-m_R';
                };
            for i=1:length(types)
                if ~nnz(strcmp(types{i},allowedTypes))
                    error(strcat(types{i},' is not a KC subtype'));
                end
            end
            indices = cell(length(types),1);
            bodyIds = cell(length(types),1);
            allKCIdsv11 = csvread('KCids_v1.1.csv');
            allNeurons = readtable('traced-neurons-v1.1.csv');
            numKCs = length(allKCIdsv11);

            for i=1:numKCs
                subtype = find(strcmp(allNeurons.instance(allNeurons.bodyId==allKCIdsv11(i)), types));
                if ~isempty(subtype)
                    bodyIds{subtype} = [bodyIds{subtype};  allKCIdsv11(i)];
                    indices{subtype} = [indices{subtype}; i-1]; %i-1 because Matlab indexes from 1 and Python indexes from 0
                end
            end
            
        end
        
        function kc = createKC(i)
            % read in the KC's skeleton
            kc = KCskel(strcat('KCskel',num2str(i),'.csv'));
            % read in the input and output synapses of this KC and map them onto the
            % skeleton
            kc = kc.getSynSets({strcat('KCoutputs',num2str(i),'.csv'), ...
                strcat('KCinputs',num2str(i),'.csv')}, [1 2]);
            % connect the real synapse locations to where they are mapped on the
            % skeleton (there is no reason for this to be a separate function - just
            % never got round to integrating)
            kc = kc.mapOrigSynsToSynLocs;
            
            kc = kc.labelSynsByROI([1 2]);
            kc = kc.definePosteriorPeduncleEdge();
            kc.saveKCskel;
            
        end
    end
end