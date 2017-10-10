function [TSubs,xySubs]=findSubArrayPhaseCenters(T,xyPoints)
%%FINDSUBARRAYPHASECENTERS Given the locations of elements in a phased
%           array that has been broken into subarrays as well as a real
%           amplitude tapering applied to each element (such as a Taylor
%           tapering for the sum beam, not a [complex] difference beam
%           tapering), approximate each of the subarrays as a single
%           antenna element with a particular amplitude tapering. The
%           location of this equivalent element is the phase center of the
%           subarray. This assumes that all elements have the same
%           directional response. While such subarray approximations will
%           be reasonable in the main beam near the look direction of the
%           subarray. They will not be reasonable outside of the main beam,
%           because the approximation introduces grating lobes.
%
%INPUTS: T A numSubarraysXnumPoints matrix of positive, real amplitude
%          weights for the elements in each subarray. A weight of zero
%          means that a particular element is not in a given subarray.
% xyPoints A 2XnumPoints set of numPoints points in the aperture plane
%          corresponding to element phase center positions. For a linear
%          array, this is a 1XnumPoints set of points.
%
%OUTPUTS: TSubs A numSubarraysXnumSubarrays diagonal matrix where
%               TSubs(i,i) is the amplitude weight for the ith subarray
%               when approximating it as a single antenna element.
%        xySubs A 2XnumSubarrays (for a planar array) or a 1XnumSubarrays
%               (for a linear array) set of the phase center locations of
%               the subarrays.
%
%The algorithm simply matches the subarray output and the gradient of the
%subarray output (taken with respect to the target location) in the look
%direction. It assumes that the response of a particular element to a
%signal in direction uVec=[u;v], where u and v are direction cosines, is
%f(u)*b*exp(-1j*2*pi*r'*(uVec-u0)), where u0 is the look direction, b is a
%complex amplitude r is the location of the element from the phase center
%of the full array (i.e. a value from xyPoints), and f(u) is the antenna
%response, which is only a positive value, but can depend on u.
%
%This function implements the approximation as derived in [1].
%
%REFERENCES:
%[1] U. Nickel, "Subarray configurations for digital beamforming with
%    low sidelobes, adaptive interference suppression and superresolution,"
%    1995, FFM-Report (available from Fraunhofer FKIE, Wachtberg, Germany).
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=size(xyPoints,1);
numSubarrays=size(T,1);

xySubs=zeros(numDims,numSubarrays);

TSubs=sum(T,2);

for curSub=1:numSubarrays
    for curDim=1:numDims
        xySubs(curDim,curSub)=sum(T(curSub,:).*xyPoints(curDim,:));
    end
    xySubs(:,curSub)=xySubs(:,curSub)/TSubs(curSub);
end

TSubs=diag(TSubs);

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
