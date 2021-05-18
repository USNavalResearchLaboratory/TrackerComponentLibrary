function circVar=findCircVar(x,w)
%%FINDCIRCVAR Given a set of samples of a circular distribution, determine
%             the circular variance. As in Chapter 2.3.1 of [1], this is
%             just 1 minus the resultant length of the first trigonometric
%             moment.
%INPUTS: x The 1XN or NX1 vector of possibly weighted samples.
%        w The NX1 weights associated with the samples. If this parameter
%          is omitted or an empty matrix is passed, then the samples are
%          uniformly weighted.
%
%OUTPUTS: circVar The circular variance of the measurements.
%
%The mean resultant length is found using findTrigMomentFromSamp and is
%then subtracted from 1.
%
%REFERENCES
%[1] K. V. Mardia and P. E. Jupp, Directional Statistics. Chichester: John
%    Wiley and Sons, 2000.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    w=[];
end

[~,~,RS]=findTrigMomentFromSamp(1,x,w);
circVar=1-RS;

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
