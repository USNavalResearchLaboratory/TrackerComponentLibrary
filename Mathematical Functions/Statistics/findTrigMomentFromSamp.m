function [rho,theta,R]=findTrigMomentFromSamp(n,x,w)
%%FINDTRIGMOMENTFROMSAMPLE Given a set of samples of a circular
%          distribution, determine a particular trigonometric sample
%          moment.
%
%INPUTS: n The order of the moment desired. This is >=1.
%        x The 1XN or NX1 vector of possibly weighted samples.
%        w The NX1 weights associated with the samples. If this parameter
%          is omitted or an empty matrix is passed, then the samples are
%          uniformly weighted.
%
%OUTPUTS: rho The (complex) mean resultant value for the nth moment. Note
%             that rho=R*exp(1j*theta).
%       theta The (real) trigonometric mean angle in radians for the nth
%             moment. This is between -pi and pi.
%           R The (real) mean resultant length for the nth moment.
%
%An expression for the complex trigonometric moments is given in Equation
%14 of [1]. This is just a weighted version of the definition in Chapter
%2.4 of [2].
%
%REFERENCES:
%[1] G. Kurz, I. Gilitschenski, R. Y. Siegwart, and U. D. Hanebeck,
%    "Methods for deterministic approximation of circular densities,"
%    Journal of Advances in Information Fusion, vol. 11, no. 2, pp.
%    138-156, Dec. 2016.
%[2] K. V. Mardia and P. E. Jupp, Directional Statistics. Chichester: John
%    Wiley and Sons, 2000.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(x);
if(nargin<3||isempty(w))
    w=ones(N,1);
end

%Normalize.
w=w/sum(w);

rho=sum(exp(1i*n*x(:)).*w(:));
if(nargout>1)
    theta=angle(rho);
    R=abs(rho);
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
