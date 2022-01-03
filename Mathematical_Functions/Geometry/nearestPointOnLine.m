function [pc,t]=nearestPointOnLine(z,p0,d)
%%NEARESTPOINTONLINE Given a point z, determine the nearest point on the
%           line given parametrically by the equation zt=p0+t*d, where t is
%           the scalar parameteric parameter. All quantities are real.
%
%INPUTS: z The numDimX1 point.
%       p0 The numDimX1 origin of the parameteric line.
%        d The numDimX1 direction of the parameteric line.
%
%OUTPUTS: pc The numDimX1 point on the line that is closest to z.
%          t The value of the parameteric parameter that corresponds to the
%            point pc on the line.
%
%The expression for the nearest point on a line is given in Chapter 10.2 of
%[1].
%
%EXAMPLE:
% p0=[0;1;0];
% d=[1;5;6];
% z=[1;36;-8];
% [pc,t]=nearestPointOnLine(z,p0,d)
% %One can verify that 
% (z-pc)'*d
% %Is zero within finite precision bounds, which is consistent with the
% %line joining the point to the specified line being orthogonal to the
% %specified line.
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

u=d/(d'*d);
t=u'*(z-p0);
pc=p0+t*d;

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
