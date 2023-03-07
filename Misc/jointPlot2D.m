function figH = jointPlot2D(x,y,pdf,jptype,figH)
%%JOINTPLOT2D A function for generating joint density plots for bivariate
%             random vectors. A joint density plot is created from [x,y]
%             pairs and colored according to the pdf input. Histograms of
%             the marginal pdfs are generated for the horizontal and
%             vertical axes.
%
%INPUTS:
% x: A vector of real-values sampled from the distribution to be plotted
%    along the horizontal axis.
% y: A vector of the same size as x containing real-values sampled from
%    the distribution to be plotted along the vertical axis.
% pdf: A vector of the same size as x defining relative intensities of
%    colors to be used.
% jptype: A string indicating the type of plot to show for the joint plot.
%         An empty string equates to the default.
%         Options are:
%         -"scatter": generates a scatter plot (default)
%         -"trisurf": generates a surface plot over a Delaunay
%                     triangulation of the data
%         -"trimesh": generates a mesh plot over a Delaunay
%                     triangulation of the data
%         If the provided option is not known, an error is returned.
% figH: An optional figure handle which will be used for holding the
%       plots. If not provided, a new figure will be generated.
%
%OUTPUTS:
% figH: A handle to the generated figure.
%
%October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if nargin<4
    jptype = [];
end
if nargin<5
    figH = figure();
else
    figure(figH)
    clf
end

%Get positions for setting based on main joint density subplot
subplot(2,2,3);
posJ = get(gca,'position'); %[left, bottom, width, height]
posJ(1,3) = posJ(1,3)*1.5; %increase the width
posJ(1,4) = posJ(1,3); %increase the height
set(gca,'position',posJ)

%Get axis limits
xrg = [min(x)*1.1,max(x)*1.1];
yrg = [min(y)*1.1,max(y)*1.1];

%Horizontal axis histogram
subplot(2,2,1)
histogram(x,'Normalization','probability')
pos = get(gca,'position'); %[left, bottom, width, height]
pos(1,4) = pos(1,4)/2; %halve the height
pos(1,2) = 1.15*(posJ(1,2)+posJ(1,4)); %shift the plot up
pos(1,3) = posJ(1,3); %set width to joint density width
set(gca,'position',pos)
xlim(xrg)

%Vertical axis histogram
subplot(2,2,4)
histogram(y,'Normalization','probability')
set(gca,'view',[90,-90])
pos = get(gca,'position'); %[left, bottom, width, height]
pos(1,3) = pos(1,3)/2; %halve the width
pos(1,1) = 1.15*(posJ(1,1)+posJ(1,3)); %shift the plot right
pos(1,4) = posJ(1,4); %set height to joint density height
set(gca,'position',pos)
xlim(yrg)

%Joint density
subplot(2,2,3);
posJ = get(gca,'position'); %[left, bottom, width, height]
posJ(1,3) = posJ(1,3)*1.5; %increase the width
posJ(1,4) = posJ(1,4)*1.5; %increase the height
set(gca,'position',posJ)
if isempty(jptype)||strcmpi(jptype,"scatter")
    scatter(x,y,[],pdf)
elseif strcmpi(jptype,"trisurf")
    dt = delaunay(x',y');
    trisurf(dt,x,y,pdf)
    set(gca,'view',[0,90])
elseif strcmpi(jptype,"trimesh")
    dt = delaunay(x',y');
    trimesh(dt,x,y,pdf)
    set(gca,'view',[0,90])
else
    error('Unknown plot type provided.')
end
xlim(xrg)
ylim(yrg)
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
