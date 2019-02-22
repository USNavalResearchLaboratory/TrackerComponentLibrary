function vals=BSplineInterpVal(x,t,a,k)
%%BSPLINEINTERPVAL Given a set of b-spline interpolation weights a and
%            knots t for scalar interpolation, interpolate values at the
%            points in x. The points x must be real, but the function being
%            interpolated can be complex.
%
%INPUTS: x The real, scalar point or matrix of real points at which one
%          wants the interpolated values.
%        t The NX1 or 1XN set of knots for the interpolation function. The
%          first and last k-1 knots are outside of the ends of the region
%          with data or mark the ends of the region with data. These values
%          are real.
%        a The npXnumSets collection of b-spline coefficients for
%          interpolating over a certain region for numSets sets of
%          interpolation problems. Such coefficients could be obtained from
%          the function BSplinePolyCoeffs along with t if one wishes to
%          interpolate over gridded data. These values can be complex.
%        k The order of the B-spline polynomials. The order of the
%          interpolation accuracy is k-1. This value is real.
%
%OUTPUTS: vals The numSetsXnumPoints matrix of interpolated values. This
%              If a is complex, then vals will be complex.
%
%Equation 7 in Chapter IX of [1] expresses an interpolated function in
%terms of a weigted combination of b-spline polynomials. This function
%implements that equation, using evalBSplinePolys to obtain sets of the
%nonzero b-spline polynomials.
%
%Examples of using this function are given in the comments to the function
%BSplinePolyFit.
%
%The first and last elements in t are not used in this function. However,
%the equations in [1] assume that those elements are present as they end up
%being used when computing integrals.
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(t);
numSets=size(a,2);

numCentPoints=N-2*(k-1)-1;
tCentral=t(k:(k+numCentPoints-1));

numVals=numel(x);
vals=zeros(numSets,numVals);
for curVal=1:numVals
    xCur=x(curVal);

    %We have to find the lowest index in tCentral  such that
    %t(idx)<=xCur<=t(idx+1). However, if we hit a repeated knot, we take
    %the last repeated value.
    [~,idx]=binSearch(tCentral,xCur,1);
    while(idx~=numCentPoints&&tCentral(idx)==tCentral(idx+1))
        idx=idx+1; 
    end

    B=evalBSplinePolys(t,xCur,k,idx+k-1);

    vals(:,curVal)=B'*a(idx:(idx+k-1),:);
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
