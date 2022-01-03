classdef DynamicPlot < handle
%%DYNAMICPLOT A class for maintaining a group of related axis tiles in one
%             figure so that added graphics objects are dropped after a
%             given number of updates. The tiles are indexed using a column
%             major scheme.
%
%October 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        %hist: A m-by-n cell array where (m,n) indicates the (row,column)
        % position of the axis tile in the figure. Each cell contains a
        % vector of graphics objects which are the objects created and
        % deleted within the class. Note that the plots are also editable
        % outside of the class, so the objects in hist are not necessarily
        % the only objects in a plot. To find objects which are currently
        % plotted, look in the corresponding histCount arrays for entries
        % which are not -Inf. Deleted objects will persist in the array
        % until overwritten to avoid unnecessary time clearing each deleted
        % object.
        hist
        
        %fig: A figure handle for the figure maintained by the class.
        fig
        
        %histLength: A nonnegative integer which determines the number of
        % calls to the advancedHist function which can be performed before
        % the object is deleted after its initial creation.
        histLength
        
        %histCount: A m-by-n cell array where (m,n) indicates the
        % (row,column) position of the axis tile in the figure. Each cell
        % contains a vector of integers indicating how many calls to
        % advanceHist remain before the object is deleted. A value of -Inf
        % indicates the corresponding object in hist has been deleted.
        histCount
        
        % m: A positive integer indicating the number of rows in the grid.
        m
        
        % n: A positive integer indicating the number of columns in the
        n
    end
    
    methods
        function this = DynamicPlot(m,n,xLab,yLab,axLims,axTitles,histLength,numObjs)
        %%DYNAMICPLOT Constructor which initializes a figure window with a
        %             grid of m-by-n tiles. Each tile is labeled using the
        %             inputs and a history length is set.
        %
        %INPUT:
        % m: A positive integer indicating the number of rows in the grid.
        % n: A positive integer indicating the number of columns in the
        %    grid.
        % xLab: A string for labeling the x-axis of every tile.
        % yLab: A string for labeling the y-axis of every tile.
        % axLims: A 1-by-4 vector of the form [xmin xmax ymin ymax] for
        %         setting the axis limits on each tile.
        % histLength: A non-negative integer indicating the number of times
        %             the advanceHist function should be called before a
        %             managed graphics object is deleted.
        % numObjs: A finite, positive integer specifying how many graphics
        %          objects to maintain. This is used to initialize the hist
        %          and histCount arrays, avoiding continual resizing of the
        %          arrays during runtime which can lead to large time
        %          costs.
        %
        %October 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            this.fig = figure();
            tiledlayout(m,n,'TileIndexing','columnmajor');
            for p = 1:m*n
                nexttile()
                xlabel(xLab)
                ylabel(yLab)
                axis(axLims)
                title(axTitles(p))
            end
            
            this.hist = cell(m,n);
            this.histCount = cell(m,n);
            for row = 1:m
                for col = 1:n
                    this.hist{row,col} = gobjects(1,numObjs);
                    this.histCount{row,col} = -Inf([1,numObjs]);
                end
            end
            
            this.histLength = histLength;
            this.m = m;
            this.n = n;
        end
        
        function addLine(this,p,lineParams)
        %%ADDLINE Adds a line object to the list of managed objects for
        %         tile p.
        %
        %INPUT:
        % p: A number indicating which tile's list to add the line object
        %    to.
        % lineParams: A cell array containing all desired line parameters.
        %             This is unpacked into the line command:
        %             line(lineParams{:}).
        %
        %October 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            if ~isempty(lineParams{1}) || ~isempty(lineParams{2})
                figure(this.fig);hold on;
                nexttile(p)
                newL = line(lineParams{:});
                [row,col] = ind2sub([this.m,this.n],p);
                idx = find(this.histCount{row,col}<0,1);
                % Check if there is an available entry. If not, replace
                % oldest.
                if isempty(idx)
                    idx = find(this.histCount{row,col}==min(this.histCount{row,col}),1);
                    delete(this.hist{row,col}(idx))
                end
                this.hist{row,col}(idx) = newL;
                this.histCount{row,col}(idx) = this.histLength;
            end
        end
        
        function advanceHist(this)
        %%ADVANCEHIST Primary function for updating the managed figure.
        %             This function decrements the histCount entries by 1,
        %             then deletes any objects with a negative count.
        %
        %October 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            shape = size(this.hist);
            for row = 1:shape(1)
                for col = 1:shape(2)
                    this.histCount{row,col} = this.histCount{row,col}-1;
                    garr = this.hist{row,col}(this.histCount{row,col}<0);
                    delete(garr);
                    this.histCount{row,col}(this.histCount{row,col}<0) = -Inf;
                end
            end
        end
    end
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.