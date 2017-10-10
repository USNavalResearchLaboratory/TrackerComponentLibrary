function vals=BSplineInterpValMultiDim(x,t,tLength,a,k)
%%BSPLINEINTERPVAL Given a hypermatrix of multivariate b-spline
%            interpolation weights a and the collection of knot components
%            t, interpolate values at the points in x.
%
%INPUTS: x The numDimXnumPoints set of points where one wishes to perform
%          interpolation. These values must be real.
%        t The maxNumKnotsXnumDims set of knots fo r the interpolation
%          function. These are actually the values in each coordinate
%          dimension given in ascending order. The full knots are implied
%          [t1,t2,etc.]=ndgrid(t(1:tLengths(1),1),t(1:tLengths(2),2),...)
%          and the dimensions can be put together into the full set of
%          numDimsXtotalNumPoints points as tTotal=[t1(:).';t2(:).';...].
%          The first and last k(curDim)-1 knots in each dimension are
%          outside of the ends of the region with data or mark the ends of
%          the region with data. The number of rows is the maximum number
%          needed for all of the dimensions. tLength says how many items
%          in each column are actually used. These values must be real.
%  tLength A numDimX1 vector where tLength(i) says the number of elements
%          in t(:,i).
%        a A hypermatrix with numDim+1 indices containing the set of
%          coefficients for the b-splines covering all of the interpolation
%          intervals for all of the sets of interpolation values. If there
%          is only one set, the final dimension is unitary meaning that
%          there are effectively only numDim dimensions. These values can
%          be real or complex.
%        k A numDimsX1 or 1XnumDims set of the order of the b-splines in
%          each dimension. The value k-1 is the polynomial order of the
%          approximation being performed. If the order is the same in all
%          dimensions, then a single scalar can be passed.
%
%OUTPUTS: vals The numSetsXnumPoints sets of interpolated points for each
%              of the numSets of interpolation coefficients given in a.
%
%Equation 7 in Chapter IX of [1] expresses an interpolated function in
%terms of a weigted combination of b-spline polynomials. Chapter 17
%discusses the multivariate extension using tensor product splines, which
%is implemented here.
%
%Examples of using this function are given in the comments to the function 
%BSplinePolyFitMultiDim.
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=size(x,1);

numA=size(a);
if(length(numA)==numDims)
    numSets=1;
else
    numSets=numA(end);
end

if(isscalar(k))
    k=repmat(k,[numDims,1]); 
end
k=k(:);

%This is the number of central points in each dimension. In one dimension,
%we would use these to select the central region of t as
% tCentral=t(k:(k+numCentPoints-1));
numCentPoints=tLength-2*(k-1)-1;
numCentPoints=numCentPoints(:);

numVals=size(x,2);
vals=zeros(numSets,numVals);
aSel=cell(1,numDims+1);%Allocate outside the loop
idx=zeros(numDims,1);%Allocate outside the loop
for curVal=1:numVals
    xCur=x(:,curVal);
    
    %In each dimension, we must find the lowest index in the central
    %region of t for that dimension such that t(idx)<=xCur<=t(idx+1).
    %However, if we hit a repeated knot, we take the last repeated value.
    
    for curDim=1:numDims
        tCentral=t(k(curDim):(k(curDim)+numCentPoints(curDim)-1),curDim);
        [~,curIdx]=binSearch(tCentral,xCur(curDim),1);
        while(curIdx~=numCentPoints(curDim)&&tCentral(curIdx)==tCentral(curIdx+1))
            curIdx=curIdx+1; 
        end
        idx(curDim)=curIdx;
    end
    
    B=evalBSplinePolysMultiDim(t,tLength,xCur,k,idx+k-1);

    for curSet=1:numSets
        for curDim=1:numDims
            aSel{curDim}=idx(curDim):(idx(curDim)+k(curDim)-1);
        end
        aSel{numDims+1}=curSet;

        temp=B.*a(aSel{:});
        vals(curSet,curVal)=sum(temp(:));
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
