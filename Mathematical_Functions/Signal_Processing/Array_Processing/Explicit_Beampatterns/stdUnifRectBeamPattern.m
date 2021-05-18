function B=stdUnifRectBeamPattern(uv,numRows,numCols)
%%STDUNIFRECTBEAMPATTERN Obtain the value of the beampattern for a
%       rectangular array consisting of uniformly weighted isotropic
%       elements arranged in a rectangular grid and spaced a half-
%       wavelength apart. As defined in Chapter 2.2 of [1], the beampattern
%       is the frequency-wavenumber response function evaluated versus
%       direction. It describes the complex gain of the array to an input
%       plane wave. The values returned by this function are real.
%
%INPUTS: uv A 2XN set of N directions offset from the pointing direction of
%           the array in terms of direction cosines. sum(uv.^2,1) should
%           all be less than or equal to 1 for results to be valid. No
%           warnings or errors are given for values outside this range.
%   numRows The number of rows of element in the array.
%   numCols The number of columns of elements in the array.
%
%OUTPUTS: B The beampattern evaluated at the points in uv.
%
%The beam pattern of a rectangular array of uniformly weighted isotropic
%antennas is given in terms of sine functions in Equation 4.26 of [1]. The
%beampattern plays a roal in determining the gain in the radar range
%equation, which itself plays a role in determining the detection
%probability of a target.
%
%EXAMPLE:
%Here, we will plot the beam pattern magnitude power in decibels over the
%viewable region.
% numRows=40;
% numCols=20;
% numPoints=100;
% uVals=linspace(-1,1,numPoints);
% [u,v]=meshgrid(uVals,uVals);
% uv=[u(:).';v(:).'];
% B=stdUnifRectBeamPattern(uv,numRows,numCols);
% %Set the values of the response outside of the valid region to zero.
% sel=sum(uv.^2,1)>1;
% B(sel)=0;
% B=reshape(B,numPoints,numPoints);
% 
% figure(1)
% clf
% surf(u,v,20*log10(abs(B)),'EdgeColor','None')
% axis([-1 1 -1 1 -50 0])
% caxis([-50 0])
% colormap(jet(256))
% view(20,45)
% light()
% h1=xlabel('u');
% h2=ylabel('v');
% h3=zlabel('Array Response');
% title('Array Power Response in Decibels')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] H. L. Van Trees, Optimum Array Processing. New York:
%    Wiley-Interscience, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. 

M=numRows;
N=numCols;

rowVals=(1/M)*sin(pi*(M/2)*uv(1,:))./sin((pi/2)*uv(1,:));
rowVals(~isfinite(rowVals))=1;%Remove singularity

colVals=(1/N)*sin(pi*(N/2)*uv(2,:))./sin((pi/2)*uv(2,:));
colVals(~isfinite(colVals))=1;%Remove singularity

B=rowVals.*colVals;
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
