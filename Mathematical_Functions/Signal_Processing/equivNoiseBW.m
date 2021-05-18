function ENBW=equivNoiseBW(w)
%%EQUIVNOISEBW Given a one-dimensional uniformely sampled window (tapering)
%              function, compute the equivalent noise bandwidth. This is
%              the width of a rectangular filter with the same peak power
%              gain having the same noise power. The equivalent noise
%              bandwidth is defined in Section IVA of [1].
%
%INPUTS: w An NX1 or 1XN all positive real window function.
%
%OUTPUTS: ENBW The value of the equivalent noise bandwidth normalized by
%              the noise power per frequency bin.
%
%The expression is Equation 11, but it should actually be multiplied by the
%length of the filter, which is done here.
%
%EXAMPLE 1:
% ENBW=equivNoiseBW(windowFunSym(1000,'triangular'))
%One will get a value of about 1.3333.
%
%EXAMPLE 2:
% ENBW=equivNoiseBW(windowFunSym(1000,'rectangular'))
%Onw will get a value of 1.
%
%REFERENCES:
%[1] F. J. Harris, "On the use of windows for harmonic analysis with the
%    discrete Fourier transform," Proceedings of the IEEE, vol. 66, no. 1,
%    pp. 172-204, Jan. 1978.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

ENBW=length(w)*sum(w.^2)/(sum(w).^2);

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
