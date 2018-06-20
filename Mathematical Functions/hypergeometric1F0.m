function S=hypergeometric1F0(a,z)
%%HYPERGEOMETRIC1F0 Evaluate the confluent hypergeometric function 1F0 for
%                   real or complex parameters.
%
%INPUTS: a A real or complex scalar value.
%        z A real or complex scalar value.
%
%OUTPUTS: S The value 1F0(a;;z).
%
%As in Equation 2.3 of [1], the hypergeometric function 1F0 is just
%(1-z)^(-a).
%
%REFERENCES:
%[1] J. Pearson, "Computation of hypergeometric functions," Master's
%    thesis, University of Oxford, Worcester College, 4 Sep. 2009.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

S=(1-z)^(-a);

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
