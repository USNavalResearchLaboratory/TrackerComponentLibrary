function [z,ze]=doubleLengthMult(x,xe,y,ye)
%%DOUBLELENGTHMULT Given two doublelength values, which could, for example,
%       be returned by exactpairMult or exactPairSum, find their product as
%       another doublelength value.
%
%INPUTS: x,xe A real doublelength floating point number. The exact number
%              is x+xe (added with infinite precision) and
%              abs(xe)<=abs(x+xe)*2^(-t)/(1+2^(-t)), where t is the number
%              of bits in the mantissa, which is 53, since it is assumed
%              that these values are floating point doubles.
%        y, ye A second real doublelength floating point number that should
%              be multiplied by to x,xe.
%
%OUTPUTS: z, ze The doublelength floating point number that is the product
%               of (x,xe) and (y,ye). If the inputs satisfy the proper
%               splitting of the number into two parts, (the relation in
%               size between abs(xe) and abs(x+xe) and similarly for y and
%               ye), then this will also satisfy the splitting of the
%               numbers. This does not necessarily hold if denormalized
%               numbers are encountered.
%
%This function implements the Mul2 algorithm of [1].
%
%Note that this assumes that the processor rounding mode has been set to
%round to "nearest," which is the default in Matlab, but which can be
%changed using the function setProcRoundingMode.
%
%REFERENCES:
%[1] T. J. Dekker, "A Floating Point Technique for Extending the Available
%    Precision," Numerische Mathematik, vol. 18, no. 3, Jun. 1971, pp.
%    224-242.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[c,ce]=exactPairMult(x,y);

ce=x*ye+xe*y+ce;
z=c+ce;
ze=c-z+ce;

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
