function r=getRange(xCart,useHalfRange,zTx,zRx)
%%GETRANGE Obtain the bistatic (or monostatic) range measurements of
%          targets in the absence of refraction under non-relativistic
%          mechanics, ignoring atmospheric effects. The transmitter and the
%          target can be collocated.
%
%INPUTS: xCart The numDimXN set of N target positions. numDim is typically
%              2 or 3.
% useHalfRange A boolean value specifying whether the bistatic range value
%             should be divided by two. This normally comes up when
%             operating in monostatic mode, so that the range reported is
%             a one-way range. The default if this parameter is not
%             provided is false.
%         zTx A numDimXN matrix of the positions of the transmitters . If
%             this parameter is omitted or an empty matrix is passed, the
%             transmitters are assumed to be at the origin. If only a
%             single vector is passed, then the transmitter position is
%             assumed the same for all of the target states being converted.
%         zRx A numDimXN matrix of the positions of the receivers. If this
%             parameter is omitted or an empty matrix is passed, the
%             receivers are assumed to be at the origin. If only a single
%             vector is passed, then the receiver position is assumed the
%             same for all of the target states being converted.
%
%OUTPUTS: r The 1XN bistatic ranges of the targets.
%
%Bistatic range is discussed in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(xCart,2);
numDim=size(xCart,1);

if(nargin<4||isempty(zRx))
    zRx=zeros(numDim,N);
elseif(size(zRx,2)==1)
    zRx=repmat(zRx,[1,N]);
end

if(nargin<3||isempty(zTx))
    zTx=zeros(numDim,N);
elseif(size(zTx,2)==1)
    zTx=repmat(zTx,[1,N]);
end

if(nargin<2||isempty(useHalfRange))
    useHalfRange=false;
end

r=sqrt(sum((xCart-zTx(1:numDim,:)).^2,1))+sqrt(sum((xCart-zRx(1:numDim,:)).^2,1));

if(useHalfRange)
   r=r/2; 
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
