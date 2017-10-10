function T=findDartboardSubarrays(aperFracs,subarraysPerRing,angShifts,xyPoints)
%%FINDDARTBOARDSUBARRAYS When designing subarrays for a circular aperture,
%          it is noted in [1] that a dartboard pattern tends to work pretty
%          well. It is also suggested that a subarray pattern returned by
%          findBaylissSubarrays be used as a rule of thumb for choosing how
%          many rings and the number of subarrays per ring in the
%          dartboard. When given the positions of elements in a 
%
%INPUTS: aperFracs A numRingsX1 or 1XnumRings vector sorted in increasing
%                  order where aperFrac(i) is the fraction of the aperture
%                  size where the ring specified ends. aperFrac(end) must
%                  equal 1 (100% of the aperture size). The aperture width
%                  is taken to be the maximum distance of a point in
%                  xyPoints from the origin.
% subarraysPerRing A numRingsX1 or 1XnumRings vector indicating the number
%                  of subarrays in each ring. These must be positive
%                  integers.
%        angShifts A numRingsX1 or 1XnumRings vector giving angular offsets
%                  in radians by which the first subarray in each ring
%                  begins. This is here to keep the rings from all lining
%                  up, which should improve how the grating lobes of the
%                  subarrays line up. If an empty matrix is passed, then
%                  angular shifts of zero are used.
%         xyPoints A 2XnumPoints set of numPoints points in the aperture
%                  plane corresponding to element positions. The points
%                  should be shifted such that the center of the aperture
%                  is the origin. The width of the aperture is the maximum
%                  distance of any point from the center.
%
%OUTPUTS: T A numSubarraysXnumPoints boolean matrix (a matrix such that
%           T(i,:) is a set of boolean values indicating which elements are
%           in each subarray). The subarrays do not overlap.
%
%EXAMPLE:
%Here, we lay out the elements for a circular array and then we break them
%into a dartboard configuration. We then display the elements with
%different symbols and colors for the subarrays.
% %First, we create a circular array. The element locations are given in
% %terms of the wavelength lambda, so lambda will not appear in the
% %equations for the sum beam.
% xyVals=getShaped2DLattice([49;49],'circular');
% 
% aperFracs=[1/5;2/5;3/5;4/5;1];
% subarraysPerRing=[4;4;8;8;8];
% degShifts=[0;2*pi/(2*4);0;2*pi/(2*8);0];
% T=findDartboardSubarrays(aperFracs,subarraysPerRing,degShifts,xyVals);
% 
% figure()
% clf
% hold on
% axis square
% numSubarrays=size(T,1);
% dispOpts={'.b','.g','.m','.c','.y','.k','.b','.r',...
%     'xm','xc','xy','xk','xb','xr','xb','xg',...
%     'ob','og','om','oc','oy','ok','ob','or'...
%     '+m','+c','+m','+y','+k','+r','+b','+g'...
%     'sb','sg','sm','sc','sm','sy','sk','sr'};
% for curSubarray=1:numSubarrays
%     points=xyVals(:,T(curSubarray,:)==1);
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
%[1] U. Nickel, "Subarray configurations for digital beamforming with
%    low sidelobes, adaptive interference suppression and superresolution,"
%    1995, FFM-Report (available from Fraunhofer FKIE, Wachtberg, Germany).
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numRings=length(aperFracs);
numSubarrays=sum(subarraysPerRing);
numEls=size(xyPoints,2);

if(aperFracs(numRings)~=1)
   error('The last element in  aperFracs must be one');
end

if(isempty(angShifts))
    angShifts=zeros(numRings,1);
end

%The maximum squared distance from the origin to a point is taken to be the
%radius of the aperture squared.
a=sqrt(max(sum(xyPoints.^2,1)));

temp=cumsum(subarraysPerRing(:));
numSubsInPriorRings=[0;temp(1:(numRings-1))];

T=false(numSubarrays,numEls);
for curEl=1:numEls
    %The fractional distance from the outer edge of the aperture where the
    %element is located.
    fracCur=norm(xyPoints(:,curEl))/a;
    
    %Determine which ring the element is in.
    [~, ringIdx]=binSearch(aperFracs,fracCur,2);
    
    %Each ring can have multiple subarrays. The subarrays are uniformly
    %spaced in angle, but the starting angle is not taken to be zero but
    %rather degShifts(ringIdx) radians.
    angle=atan2(xyPoints(2,curEl),xyPoints(1,curEl));
    
    subArrayWidth=2*pi/subarraysPerRing(ringIdx);
    
    subArrayInRing=fix(wrapRange(angle+angShifts(ringIdx),0,2*pi)/subArrayWidth)+1;
    subArrayIdx=numSubsInPriorRings(ringIdx)+subArrayInRing;
    T(subArrayIdx,curEl)=1;
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
