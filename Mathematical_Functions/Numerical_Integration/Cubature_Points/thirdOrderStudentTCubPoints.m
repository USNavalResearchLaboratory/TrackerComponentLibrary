function [xi,w]=thirdOrderStudentTCubPoints(mu,SR,nu)
%%THIRDORDERSTUDENTTCUBPOINTS Generate third order cubature points for
%          integration involving a multidimensional Student-t distribution.
%
%INPUTS: mu The numDimX1 mean of the distribution.
%        SR The lower-triangular square root of the numDimXnumDim positive
%           definite scale matrix of the distribution.
%        nu The scalar number of degrees of freedom of the distribution
%           nu>2.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The algorithm implemented is taken from [1]. To use Student-t measurement
%or process noise in a cubature filtering algorithm, simply replace the
%Gaussian cubature points with Student-t cubature points.
%
%REFERENCES:
%[1] Y. Huang, Y. Zhang, N. Li, S. M. Naqvi, and J. Chambers, "A robust
%    Student's t based cubature filter," in Proceedings of the 19th
%    International Conference on Information Fusion, Heidelberg, Germany,
%    5-8 Jul. 2016.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(mu,1);

w=(1/(2*n))*ones(2*n,1);
xi=zeros(n,2*n);

nuConst=sqrt(nu*n/(nu-2));

curPoint=1;
for curDim=1:n
    xi(:,curPoint)=nuConst*SR(:,curDim);
    curPoint=curPoint+1;
    xi(:,curPoint)=-xi(:,curPoint-1);
    curPoint=curPoint+1;
end

xi=bsxfun(@plus,mu,xi);
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
