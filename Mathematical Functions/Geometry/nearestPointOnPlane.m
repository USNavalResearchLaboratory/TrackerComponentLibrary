function zp=nearestPointOnPlane(z,n,p0)
%%NEARESTPOINTONPLANE Given a point z, determine the nearest point on the
%           plane given by the equation n'*(x-p0)=0, where n is a normal to
%           the plane, p0 is a point and x is the free variable. All
%           quantities are real.
%
%INPUTS: z A numDimXnumPoints set of points to project.
%        n A numDimX1 normal to the plane.
%       p0 A numDimX1 point on the plane.
%
%OUTPUTS: zp The numDimXnumPoints set of projections of the points in z
%            onto the plane.
%
%The expression for the nearest point to a plane is given in Chapter 12.1
%of [1].
%
%EXAMPLE:
% n=[1;2;3];
% p0=[0;-1;2];
% z=[4;4;4];
% zp=nearestPointOnPlane(z,n,p0)
% %One will get zp=[2.571428571428572;1.142857142857143;-0.285714285714286].
% %It can be verified that 
% n'*(zp-p0)
% %is zero within finite precision limits and 
% cross(n,(z-zp))
% %is also zero within finite precision limits.
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

d=-n'*p0;
%The equation for the plane is now:
%x'*n+d=0;

u=n/(n'*n);

%For a single point z, the following line is equivalent to
%zp=z-(z'*n+d)*u;
zp=bsxfun(@minus,z,bsxfun(@times,sum(bsxfun(@times,z,n),1)+d,u));

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
