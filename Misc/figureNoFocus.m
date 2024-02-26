function h=figureNoFocus(h)
%%FIGURENOFOCUS Select a figure to which plot commands will go, but do not
% bring the figure wuindow to the front. A figure handle or window number
% can be passed. If nothing or a number of a non-existent window is passed,
% then the window is created, which will bring it to the front. When
% plotting across multiple windows, using this function instead of the
% figure function helps prevent annoying flickering of the windows back and
% forth.
% 
%INPUTS: h A Figure handle object or a figure number. Omitting this term or
%          passing an empty matrix nothing just opens a new figure window.
%
%OUTPUTS: h This is the Figure handle object for the selected figure. If a
%           number is passed, this is the object corresponding to that
%           number.
%
%EXAMPLE:
%Here, we create windows 1 and 2 and randomly switch between plotting a
%point to window 1 and a point to window 2. Using this function in the loop
%instead of figure, we do not constantly have the windows swapping which is
%in front and the plots are made significantly faster.
% figure(1);
% clf
% hold on
% figure(2);
% clf
% hold on
% for k=1:100
%     idxSel=mod(k,2)+1;
%     figureNoFocus(idxSel);
%     z=rand(2,1)+idxSel;
%     if(idxSel==1)
%         scatter(z(1),z(2),'.b')
%     else
%         scatter(z(1),z(2),'.r')
%     end
% end
%
%August 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin>0&&~isempty(h))
    %If an input was passed, it could be a figure window number or a figure
    %window handle.
    if(ishandle(h))
        %If we are here, h is either a valid figure window handle, or the
        %valid number of a window. We wish to select the window without

        %Plot commands will not go to this figure.
        set(0, 'CurrentFigure', h); % remove focus but addressing this figure with plot commands

        %Get a figure handle to return if h was an integer.
        if (isnumeric(h))
            h=findobj('Type', 'figure','Number',h);
        end
    elseif(isnumeric(h)&&isreal(h)&&(h==fix(h))&&(h>0)&&h<(intmax-1))
        %A valid figure window number was passed, but no window with that
        %number exists, so this created ones, which bring the window to the
        %front.
        h=figure(h);
    else
        %The window number limit of intmax-1 is correct as of Matlab
        %2022b.
        error('The parameter h passed is invalid. 0<h<=intmax-1.')  
    end
else
    %If no inputs were passed, then create a new figure window and
    %inevitably bring the window to the front.
    h=figure();
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
