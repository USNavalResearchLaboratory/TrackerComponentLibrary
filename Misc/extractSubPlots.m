function extractSubPlots(theFig)
%%EXTRACTSUBFIGURES Given a figure that has been broken into subplots, this
%   function extracts each of the subflots into a new figure.
%
%INPUTS: theFig This can be either the scalar number of the figure, or the
%               figure handle object.
%
%OUTPUTS: None; new plots are made, one for each subfigure in the original
%         figure.
%
%EXAMPLE:
%A figure with 3 subplots is made. This function is then called to extract
%all 3 subplots into new figures.
% x=linspace(0,10,1000);
% y1=x.^2;
% y2=x+x.^3+2;
% y3=x.^4-2*x+x.^3-1;
% 
% figHandle=figure(1);
% clf
% hold on
% subplot(3,1,1)
% plot(x,y1)
% xlabel('Dogs')
% ylabel('Lice')
% subplot(3,1,2)
% [X,Y]=meshgrid(1:0.5:10,1:20);
% Z=sin(X)+cos(Y);
% surf(X,Y,Z)
% xlabel('Cookies')
% ylabel('Cake')
% zlabel('Muskrats')
% subplot(3,1,3)
% hold on
% plot(x,y2)
% plot(x,y3)
% xlabel('Rabbits')
% ylabel('Foxes')
% legend('Light','Dark','location','southwest')
% extractSubPlots(figHandle)
%
%May 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(theFig))
    error('The figure must be specified.')
end

if(isscalar(theFig))
    %If a figure number is passed.
    figHandle=figure(theFig);
elseif(isa(theFig,'matlab.ui.Figure'))
    figHandle=theFig;
else
    error('TheFig must be either a figure number or a figure handle object.')
end

theChildren=get(figHandle,'Children');
numChildren=size(theChildren,1);
childIdx=1;
while(childIdx<=numChildren)
    newFig=figure();
    curChild=theChildren(childIdx);
    if(isempty(curChild.Tag))%There is no legend.
        copyobj(curChild,newFig)
        set(newFig.Children(1),'Position', get(0, 'DefaultAxesPosition'));
        childIdx=childIdx+1;
    else
        %There is a legend, so two objects must be copied.
        copyobj([curChild,theChildren(childIdx+1)],newFig);
        set(newFig.Children(2),'Position', get(0, 'DefaultAxesPosition'));
        childIdx=childIdx+2;
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
