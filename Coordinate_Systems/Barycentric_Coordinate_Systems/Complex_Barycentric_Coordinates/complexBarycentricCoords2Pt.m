function x=complexBarycentricCoords2Pt(c,v)
%%COMPLEXBARYCENTRICCOORDS2PT Given a set of complex barycentric
%   coordinates and a st of vertices, get the real 2D point that is
%   specified by the coordinates.
%
%INPUTS: c A numVertXnumPts set of complex coordinate vectors to convert.
%        v This is either a 2XN set of real vertices, or a length N vctor
%          of complex scalar vertices (x coordinate is real and y
%          coordinate is imaginery). These define the polygon over which
%          the vertices are specified.
%
%OUTPUTS: x The 2XnumPts set of the points given by the complex barycentric
%           coordinates.
%
%Complex barycentric coordinates are discussed in [1].
%
%REFERENCES:
%[1] O. Weber, M. Ben-Chen, C. Gotsman, and K. Hormann, "A complex view of
%    barycentric mappings," Computer Graphics Forum, vol. 30, no. 5, pp.
%    1533-1542, Aug. 2011.
%
%October 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPts=size(c,2);

if(isvector(v))
    zj=v;
else
    zj=v(1,:).'+1j*v(2,:).';
end

x=zeros(2,numPts);
for curPt=1:numPts
    xComplex=sum(c(:,curPt).*zj);
    x(:,curPt)=[real(xComplex);imag(xComplex)];
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
