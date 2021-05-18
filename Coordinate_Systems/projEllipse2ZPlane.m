function [M,kpp,c]=projEllipse2ZPlane(M,kpp,c,zHeight)
%%PROJELLIPSE2ZPLANE Given a 3D ellipsoid specified in the form 
%            (t-kpp)'*M*(t-kpp)=c, where t is a 3X1 point on the surface of
%            the ellipsoid, determine the ellipse that is obtained by
%            clutting the ellipsoid by an x-y plane located with a z height
%            of zHeight. This function can be useful for making 2D plots of
%            3D ellipsoids.
%
%INPUTS: M A 3X3 symmetric matrix.
%      kpp A 3X1 vector.
%        c A non-negative scalar value.
%
%OUTPUTS: M A 2X2 symmetric matrix, or an empty matrix if the x-y plane at
%           the specified height does not intersect the ellipsoid.
%       kpp A 2X1 vector, or an empty matrix if there is no intersection.
%         c A non-negative scalar value, or an empty matrix if there is no
%           intersection.
%
%Given the initial equation (t-kpp)'*M*(t-kpp)=c, one can multiple out the
%known component in t to get a general quadratic equation that satisfies
%the input requirements to the function quadEq2CenteredForm. Thus, the
%function quadEq2CenteredForm is used to put the result into the centered
%quadratic form used for the output. If the c value returned by
%quadEq2CenteredForm is negative, that indicates that there is no
%intersection of the plane and the ellipse, so empty matrices can be
%returned. 
%
%November 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

deltaVal=zHeight-kpp(3);

A=M(1:2,1:2);
bVec=zeros(2,1);
bVec(1)=deltaVal*(M(1,3)+M(3,1));
bVec(2)=deltaVal*(M(2,3)+M(3,2));
cVal=M(3,3)*deltaVal^2-c;

[kVec,cNew]=quadEq2CenteredForm(A,bVec,cVal);

%If the plane intersects the ellipsoid.
if(cNew>=0)
    M=A;
    kpp=kVec+kpp(1:2);
    c=cNew;
else%The plane does not intersect the ellipsoid.
    M=[];
    kpp=[];
    c=[];
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
