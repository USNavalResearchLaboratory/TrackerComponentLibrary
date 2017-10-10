function vals=invGudermannian(x)
%%INVGUDERMANNIAN  Evaluate the inverse Gudermannian function of x. This is
%               the INVERSE of integral from 0 to x of (1/cosh(t)) dt. The
%               inverse Gudermannian function of a latitude on a sphere
%               provides the isometric latitude, which plays a role in the
%               transverse Mercator projection.
%
%INPUTS: x A vector or matrix of values at which the inverse Gudermannian
%          should be evaluated. x can be between -pi/2 and pi/2.
%
%OUTPUTS: val The inverse Gudermannian evaluated at each point in x. It
%             ranges from -Inf to Inf.
%
%A noted at [1], the Gudermannian can be easily expressed in terms of
%hyperbolic and regular trigonometric functions. The inverse function comes
%from inverting those.
%
%This is the inverse of the function Gudermannian.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gudermannian." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Gudermannian.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    vals=atanh(sin(x));
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
