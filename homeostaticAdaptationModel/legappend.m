function [legend_h,object_h,plot_h,text_strings] = legappend(newStrings,varargin)
%LEGAPPEND appends new entries to the end of a legend by deleting the
%current legend and recreating a new, similar legend. 
% 
%% Syntax
% 
% legappend('new legend entry') 
% legappend('new entry 1','new entry 2',...,'new entry N') 
% legappend('') 
% legappend('','',...,'')
% [legend_h,object_h,plot_h,text_strings] = legappend(...)
% 
%% Description 
% 
% legappend('new legend entry') appends an existing legend with "new
% legend entry".
% 
% legappend('new entry 1','new entry 2',...,'new entry N') adds several
% new entries to the legend. 
% 
% legappend('') deletes the last entry from the legend. 
% 
% legappend('','',...,'') deletes the last several entries from the
% legend.
% 
% [legend_h,object_h,plot_h,text_strings] = legappend(...) returns legend_h, the
% handle of the new legend; object_h, handles of the line, patch, and
% text graphics objects used in the legend; plot_h, handles of the lines
% and other objects used in the plot; and text_strings, a cell array of
% the text strings used in the legend. Note that for new legend entries,
% legappend does not add entries to a current legend, but deletes the
% current legend and recreates a new one. As a result, the legend handle
% will change with each new-entry use of legappend.  The legend handle
% does not change when legappend is used to delete an entry. 
% 
% 
%% Author Info
% Version 2 of this script was written by Samantha Kulpinski of General
% Electric, October 2020.  Example scripts were written for the original
% version of this code by Chad A. Greene.
%
% The original version of this function was created by Chad A. Greene of 
% the Institute for Geophysics, The University of Texas at Austin, 
% July 2014. 
% 
% See also legend.

v=version;
v=split(v,' ');
v=split(v{1},'.');
a = gca;

% in R2016a (version 9.0) it appears that legend changed from being a child
% of an Axes to being a child of a Figure, thus the identification of a
% legend object varies based on MATLAB version:
if v{1} >=9  % version R2016a or later
    h = findobj(gcf,'Type','Legend');

    length_h=length(h);
    if length_h == 1
        % check if we are looking at correct axis
        if ~(h.Axes == a)
            h=gobjects(0,0); % create graphics placeholder
        end
    elseif length_h > 1
        for ii=1:length_h
            if h(ii).Axes == a
                h = h(ii);
                break
            end
        end
    end
    % if no legend exists, create placeholders
    if isempty(h)
        str = {}; 
    else
        str = get(h,'String'); 
    end


    if exist('varargin','var') 
        newStrings = [newStrings,varargin];
    end
    deleteEntries = sum(cellfun('isempty',newStrings));
    if isempty(newStrings) 
        deleteEntries = 1; 
    end

    if ~deleteEntries
        if iscell(newStrings)
            for k = 1:length(newStrings) 
                str{end+1}=newStrings{k}; 
            end
        end 
        if ~iscell(newStrings)
            str{end+1}=newStrings; 
        end

        [legend_h,object_h,plot_h,text_strings] = legend(a,str);
        h = findobj(gcf,'Type','Legend');
        
        length_h=length(h);
        if length_h > 1
            for ii=1:length_h
                if h(ii).Axes == a
                    h = h(ii);
                    break
                end
            end
        end


        prop.boxon = get(h,'visible');
        prop.loc = get(h,'location'); 
        prop.color = get(h,'color'); 
        prop.orient = get(h,'Orientation'); 

        if strcmpi({prop.boxon},'off')
            legend boxoff
        end

        set(legend_h,'location',prop.loc,'color',prop.color,'Orientation',prop.orient)


    end

    if deleteEntries
        set(h,'String',str(1:end-nargin))
        [legend_h,object_h,plot_h,text_strings] = legend;
    end


    if nargout==0
        clear legend_h object_h plot_h text_strings
    end
else % earlier than R2016a, original code from first version of legappend:
    %% Original version of script
    h =  findobj(gcf,'Type','axes','Tag','legend');

    prop.boxon = get(h,'visible');
    prop.loc = get(h,'location'); 
    prop.color = get(h,'color'); 
    prop.orient = get(h,'Orientation'); 



    allDatah = flipud(get(gca,'children')); 
    str = get(h,'String'); 

    if exist('varargin','var') 
        newStrings = [newStrings,varargin];
    end
    deleteEntries = sum(cellfun('isempty',newStrings));
    if isempty(newStrings) 
        deleteEntries = 1; 
    end

    if ~deleteEntries
        if iscell(newStrings)
            for k = 1:length(newStrings) 
                str{end+1}=newStrings{k}; 
            end
        end 
        if ~iscell(newStrings)
            str{end+1}=newStrings; 
        end


        [legend_h,object_h,plot_h,text_strings] = legend(h,allDatah,str);

        if strcmpi({prop.boxon},'off')
            legend boxoff
        end

        set(legend_h,'location',prop.loc,'color',prop.color,'Orientation',prop.orient)


    end

    if deleteEntries
        set(h,'String',str(1:end-nargin))
        [legend_h,object_h,plot_h,text_strings] = legend;
    end


    if nargout==0
        clear legend_h object_h plot_h text_strings
    end

%% 
end





