function r=findSubarrayWeights(T,d,a0)
%%FINDSUBARRAYWEIGHTS Given that the elements in an array (linear, planar)
%           have been grouped into subarrays and some type of element level
%           tapering might be applied, find the subarray-level tapering
%           that best approximates a tapering desired for alll elements in
%           an array. For example, one mgiht taper all the elements with a
%           Taylor tapering to get a good sum beam, but then want to design
%           subarray weights to be able to get a good difference beam.
%
%INPUTS: T A numSubarraysXnumEls matrix holding that breaks the elements
%          into subarrays and includes any element-level tapering.
%        d A numElsX1 vector of the desired element-level tapering. This
%          function provides a weighting for the subarrays to approximate
%          this.
%       a0 An optional numElX1 parameter. If one wishes to place a null in
%          the direction of a0 (often one might want to place a null in the
%          look direction when approximating a difference beam), then this
%          input can be provided. Otherwise, unconsitrained optimization is
%          performed. This parameter is the array manifold in the look
%          direction. For boresight, this will usually just be a numELX1
%          vector of ones. For an isotropic model, the vector is typically 
%          a0=exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,u),1)).' where
%          xyPoints is a 2XnumEl matrix of the locations of the elements in
%          the array.
%
%OUTPUTS: r The set of (possibly complex) weights for the subarrays to
%           approximate the given desired element-level tapering.
%
%This function implements the equations in [1]. The unconstrained
%optimization is just r=arg min_{r} norm(T*r-d)^2. The constrained
%optimization with a0 just adds the constraint that a0'*T*r'=0.
%
%Two approaches are given in the paper. The second approach considers
%emphasizing certain directions using a penalty function. However, when
%considering approximating the tapering for a Bayliss-weighted difference
%pattern, it was shown that this method was essentially the same as the
%simpler first method. Thus, this function just implements the first
%method.
%
%REFERENCES:
%[1] U. R. O. Nickel, "Subarray configurations for digital beamforming with
%    low sidelobes and adaptive interference suppression," in Record of the
%    IEEE International Radar Conference, Alexandria, VA, 8-11 May 1995,
%    pp. 714-719.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

T=T.';

numEls=size(T,1);

%Note that pinv(T)=inv(T'*T)*T'
TInv=pinv(T);

%If the constraint on placing a null in the look direction is imposed:
if(nargin>2&&~isempty(a0))
    r=TInv*(eye(numEls,numEls)-(a0*a0'*T*TInv)/(a0'*T*TInv*a0))*d;
else
    r=TInv*d;
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
