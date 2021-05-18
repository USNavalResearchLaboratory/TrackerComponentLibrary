function xyz=LambertProj2Cart(xy,r)
%%LAMBERTPROJ2CART Map points in a 2D Lambert projection back onto a 3D
%           sphere of radius r.
%
%INPUTS: xy The 2XnumPts set of points in the Lambert projection.
%         r The radius of the sphere from which the Lambery projection was
%           derived. If omitted or an empty matrix is passed, the default
%           of r=1 is used.
%
%OUTPUTS: xyz The 3XnumPts set of points mapped onto the sphere of radius
%             r.
%
%This implements Equation 7 in [1].
%
%REFERENCES:
%[1] D. Rosca "New uniform grids on the sphere," Astronomy & Astrophysics,
%    vol. 520, no. A63, Sep.-Oct. 2010.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(r))
    r=1;
end

numPts=size(xy,2);

r2=r^2;

xyz=zeros(3,numPts);
for curPt=1:numPts
    x=xy(1,curPt);
    y=xy(2,curPt);
    
    const=sqrt(1-(x^2+y^2)/(4*r2));
    xyz(:,curPt)=[const*x;const*y;-r+(x^2+y^2)/(2*r)];    
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
