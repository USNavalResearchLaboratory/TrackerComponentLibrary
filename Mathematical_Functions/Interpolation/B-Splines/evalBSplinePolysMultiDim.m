function B=evalBSplinePolysMultiDim(t,tLength,x,k,idx)
%%EVALBSPLINEPOLYSMULTIDIM Given a set of knots and a multidimensional
%           point, evaluate all associated order-k (k veing a vector)
%           multivariate b-spline polynomials that might be nonzero. The
%           result is just the repeated outper product of all of the
%           one-dimensional ones that one can obtain using the function
%           evalBSplinePolys in one dimension. This function cannot handle
%           values of x that are outside of the region 
%           t(idx+-1k,curDim)<=x(curDim)<=t(idx+k,curDim) in any dimension.
%
%INPUTS: t The maxNumKnotsXnumDims set of knots for the interpolation
%          function. These are actually the values in each coordinate
%          dimension given in ascending order. The full knots are implied
%          [t1,t2,etc.]=ndgrid(t(1:tLengths(1),1),t(1:tLengths(2),2),...)
%          and the dimensions can be put together into the full set of
%          numDimsXtotalNumPoints points as tTotal=[t1(:).';t2(:).';...].
%          The first and last k(curDim)-1 knots in each dimension are
%          outside of the ends of the region with data or mark the ends of
%          the region with data. The number of rows is the maximum number
%          needed for all of the dimensions. tLength says how many items
%          in each column are actually used.
%  tLength A numDimX1 vector where tLength(i) says the number of elements
%          in t(:,i).
%        x The numDimXnumPoints set of points where one wishes to evaluate
%          the multivariate b-spline weight.
%        k A numDimsX1 or 1XnumDims set of the order of the b-splines in
%          each dimension. The value k-1 is the polynomial order of the
%          approximation being performed. If the order is the same in all
%          dimensions, then a single scalar can be passed.
%      idx This input specifies where in terms of the knots the B-spline
%          polynomials are centered. This is a numDimX1 integer value
%          where the ith dimension ranges from 1 to tLength(i)-k.
%
%OUTPUTS: B In the ith dimension, the nonzero scalar B-spline values go
%           from index idx(i)-k+1 to idx(i). This is a kXkXk...
%           hypermatrix. The number of indices correspond to the number of
%           dimensions and the value in the i,,j,k etc-th index is the
%           value of the tensor product of the b-splies in the
%           i,j,k,etc.-th offset from idx-k in indexation. 
%
%B-spline polynomials are discussed in detail in Chapter IX of [1]. Tensor
%product splines are discussed in Chapter 17. This function creates the
%matrix by calling evalBSplinePolys in each dimension and combining the
%outputs with repeated calls to outerProductTensor.
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);
numPoints=size(x,2);

if(isscalar(k)&&numPoints>1)
   k=repmat(k,[numPoints,1]);
end

%Find the nonzero vectors for each of the dimensions.
B=evalBSplinePolys(t(1:tLength(1),1),x(1),k(1),idx(1));
for curDim=2:numDim
    BCur=evalBSplinePolys(t(1:tLength(curDim),curDim),x(curDim),k(curDim),idx(curDim));
    B=outerProductTensor(B,BCur);
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
