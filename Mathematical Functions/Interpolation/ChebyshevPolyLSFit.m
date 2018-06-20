function a=ChebyshevPolyLSFit(tauList,n,fList,dfList,w)
%%CHEBYSHEVPOLYLSFIT  Find the weights for a series of Chebyshev
%                     polynomials of the first kind to best fit data given
%                     at a number of discrete, scalar points, constrained
%                     such that the values at the endpoints are perfectly
%                     met. If desired, derivatives can also be added to the
%                     least squares problem, whereby the derivatives at the
%                     endpoints are perfectly matched.
%
%INPUTS: tauList An NX1 or 1XN set of points for which the interpolation
%                weights should be found. tauList should be sorted in
%                increasing order.
%              n The non-negative integer maximum order of the
%                highest-order Chebyshev polynomial in the expansion.
%          fList The function to be interpolated at the values of tauList.
%         dfList An optional set of derivatives of the function to be
%                interpolated evaluated at the points in tauList.
%              w If dfList is provided, one can optionally provide w, which
%                is a 2X1 or 1X2 weighting vector, where w(1) is the
%                weighting for the squared function errors and w(2) is the
%                weighting for the squared derivative errors. Different
%                weightings can be useful in the least squared problem,
%                because the function and its derivative will typically
%                span different domains.
%
%OUTPUTS: a A set of weights that can be used with the function
%           ChebyshevPolySynth(tau,a,tauStart,tauEnd) to interpolate a
%           point tau, where tauStart=tauList(1) and tauEnd=tauList(end)
%           and the points tau is within the starting and ening  bounds.
%
%The algorithm implemented is taken from [1]. There it is suggested that
%when performing least squares Chebyshev polynomial fitting with
%ephemerides, one use w(1)=1 and w(2)=0.16.
%
%REFERENCES:
%[1] X. X. Newhall, "Numerical representation of planetary ephemerides,"
%    Celestial Mechanics, vol. 45, no. 1-3, pp. 305-310, 1989.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

NPoints=length(fList);

tauStart=tauList(1);
tauEnd=tauList(end);

%If no derivative information is provided.
if(nargin<4||isempty(dfList))
    %Each row of T is a different point. Each column of T is a different
    %order.
    T=ChebyshevPoly(tauList,n,tauStart,tauEnd)';
    
    %There are two constraints, one for each endpoint.
    constEq=[T(1,:);T(end,:)];
    C1=[T'*T,constEq'
       constEq,zeros(2,2)];

    zMat=zeros(2,NPoints);
    zMat(1,1)=1;
    zMat(2,end)=1;

    C2=[T';
        zMat];

    aLambda=C1\C2*fList(:);
    a=aLambda(1:(n+1));
else
    if(nargin<5||isempty(w))
        w=[1;1];
    end
    
    [dT,T]=ChebyshevPolyDeriv(tauList,n,tauStart,tauEnd);
    %Each row of dT and T is a different point. Each column of is a
    %different order.
    T=T';
    dT=dT';
    
    T=[T;dT];
    %Thereare four constraints, one for each endpoing and its derivative.
    constEq=[T(1,:);T(NPoints,:);T(NPoints+1,:);T(end,:)];
    
    W=diag([ones(NPoints,1)*w(1);ones(NPoints,1)*w(2)]);
    
    C1=[T'*W*T,constEq'
       constEq,zeros(4,4)];
    
    zMat=zeros(4,NPoints*2);
    zMat(1,1)=1;
    zMat(2,NPoints)=1;
    zMat(3,NPoints+1)=1;
    zMat(4,end)=1;
    C2=[T'*W;
        zMat];

    aLambda=C1\C2*[fList(:);dfList(:)];
    a=aLambda(1:(n+1));
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
