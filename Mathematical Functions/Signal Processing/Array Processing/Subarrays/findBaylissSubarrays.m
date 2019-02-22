function T=findBaylissSubarrays(xyPoints,numLevels,sidelobedB,N,adjMat)
%%FINDBAYLISSSUBARRAYS Given a elements for a linear array or for a planar
%          array with a circular aperture, group the elements into
%          subarrays that discretize Bayliss difference patterns to a given
%          number of levels. For a planar array, both horizontal and
%          vertical difference pattern are overlapped so that the subarrays
%          can produce difference patterns in either direction. This
%          implements the basic method of determining subarrays as in [1].
%
%INPUTS: xyPoints A 2XnumPoints set of numPoints points in the aperture
%             plane corresponding to element positions. For a linear array,
%             this is a 1XnumPoints set of points. The points should be
%             shifted such that the center of the aperture is the origin.
%             The width of the aperture is the maximum distance of any
%             point from the center.
%   numLevels The number of levels to use in discretizing the Bayliss
%             pattern. If this parameter is omitted or an empty matrix is
%             passed, a default of 4 is used.
%  sidelobedB The number of decibels of the ratio of the close-in
%             sidelobe voltages to the main lobe voltage in the differnce
%             pattern. This must be a negative number. A typical value is
%             -30.
%           N The Bayliss tapering is computed using a certain number of
%             terms. Using too many terms can be undesirable as edge
%             illumination increases, as noted in [1]. If this parameter is
%             omitted or an empty matrix is passed, then the default of 17
%             is used. In [1], it is suggested that N be chosen to be
%             <2*a/lambda, where a is the radius of the aperture and
%             lambda the wavelength.
%      adjMat This parameter is only used with planar arrays and is
%             optional. This is a numPointsXnumPoints adjacency matrix such
%             that adjMat(i,j) is nonzero if element i is adjacent to
%             element j. The value of adjMat(i,i) does not matter. Only the
%             lower half of this matrix is used. If this matrix is omitted
%             or an empty matrix is passed, all elements are taken to be
%             adjacent to each other. The purpose of this matrix is to make
%             sure that no disjoint subarrays are formed. That is, a
%             subarray consists of a continuum of adjacent elements with no
%             breaks.
%
%OUTPUTS: T A numSubarraysXnumPoints boolean matrix (a matrix such that
%           T(i,:) is a set of boolean values indicating which elements are
%           in each subarray). The subarrays do not overlap.
%
%The design of subarrays is a very difficult optimization problem. This
%implements the algorithm of [1], which is based on the idea that the
%ability to form good difference beams should play a role in subarray
%design. The method presented is significatnly simpler than many of the
%genetic optimization methods in the literature.
%
%EXAMPLE 1:
%This is an example of a linear array being broken into six subarrays.
% N=17;
% sidelobedB=-30;
% Nx=41;%There are 2*Nx+1 points total.
% %Generate points symmetric about the origin.
% xPoints=(-(Nx-1)/2):1/2:((Nx-1)/2);
% numLevels=6;%Discretize into four levels.
% T=findBaylissSubarrays(xPoints,numLevels,sidelobedB,N);
% 
% %Display the elements in the array with different coloring for each
% %subarray.
% figure(1)
% clf
% hold on
% els=xPoints(T(1,:));
% numEls=length(els);
% scatter(els,zeros(1,numEls),'or')
% els=xPoints(T(2,:));
% numEls=length(els);
% scatter(els,zeros(1,numEls),'xg')
% els=xPoints(T(3,:));
% numEls=length(els);
% scatter(els,zeros(1,numEls),'+b')
% els=xPoints(T(4,:));
% numEls=length(els);
% scatter(els,zeros(1,numEls),'*k')
% els=xPoints(T(5,:));
% numEls=length(els);
% scatter(els,zeros(1,numEls),'sc')
% els=xPoints(T(6,:));
% numEls=length(els);
% scatter(els,zeros(1,numEls),'dm')
%
%EXAMPLE 2:
%This is an example of a circular planar array. This requires the formation
%of two Bayliss patterns that are them discretized and overlapped. The
%elements in one row are offset by lambda/4 from those in the next, so they
%do not form a very regular grid.
% %First, we create a circular array. The element locations are given in
% %terms of the wavelength lambda, so lambda will not appear in the
% %equations for the sum beam.
% [xyVals,adjMat]=getShaped2DLattice([49;49],'circular');
% 
% sidelobedB=-30;
% N=5;
% numLevels=3;
% T=findBaylissSubarrays(xyVals,numLevels,sidelobedB,N,adjMat);
% 
% %Display the subarrays with different colors and symbols.
% figure(2)
% clf
% hold on
% axis square
% numSubarrays=size(T,1);
% dispOpts={'.b','og','xr','+c','*m','sy','dk','vg','^b','<r','>c','pm','hy'};
% for curSubarray=1:numSubarrays
%     points=xyVals(:,T(curSubarray,:));
%     optVal=mod(curSubarray,length(dispOpts))+1;
%     scatter(points(1,:),points(2,:),dispOpts{optVal},'linewidth',2);
% end
% h1=xlabel('x');
% h2=ylabel('y');
% title('Coloring Represents Subarrays')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] U. R. O. Nickel, "Subarray configurations for digital beamforming with
%    low sidelobes and adaptive interference suppression," in Record of the
%    IEEE International Radar Conference, Alexandria, VA, 8-11 May 1995,
%    pp. 714-719.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(N))
    N=17; 
end

if(nargin<3||isempty(sidelobedB))
    sidelobedB=-30; 
end

if(nargin<2||isempty(numLevels))
   numLevels=4; 
end

switch(size(xyPoints,1))
    case 1%Linear array
        %Bayliss weights are all imaginary.
        g=BaylissLinearTapering(sidelobedB,N,xyPoints);
        T=discretizeAndOverlapVals(imag(g),[],numLevels);
    case 2%Planar array
        numEls=size(xyPoints,2);
        
        if(nargin<5||isempty(adjMat))
            %If no adjacency matrix is given, then everything is adjacent.
            adjMat=ones(numEls,numEls);
        end
        
        %Bayliss weights are all complex.
        %Weights for horizontal differencing (in u)
        gHoriz=BaylissTapering(sidelobedB,N,xyPoints);
        %Weights for vertical differencing (in v)
        gVert=BaylissTapering(sidelobedB,N,flipud(xyPoints));
        
        %Bayliss weights are all imaginary.
        T=discretizeAndOverlapVals(imag(gHoriz),imag(gVert),numLevels,adjMat);
    otherwise
        error('Only linear and planar arrays are implemented.')
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
