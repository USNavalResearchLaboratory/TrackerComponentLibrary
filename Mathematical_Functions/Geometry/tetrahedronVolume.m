function [V,dVdv]=tetrahedronVolume(v,onlyPositive,derivIdx)
%%TETRAHEDRONVOLUME Find the signed volume volume of a tetrahedron in 3D
%                   given 4 points. The sign of the volume depends on the
%                   ordering of the points.
%
%INPUTS: v A 3X4 set of 4 vertices of the tetrahedron.
% onlyPositive If one wishes for the volume to always be positive, then
%          this value should be true. The default if omitted or an empty
%          matrix is passed is false.
% derivIdx If this function returns derivatives, this is the index of the
%          vector in v with respect to which the derivatives are taken.
%
%OUTPUTS: V The signed volume of the tetrahedron (or positive volume if
%           onlyPositive is true).
%     dVdv The 1X3 derivative of the positive area with respect to the
%           elements of v(:,derivIdx) in the order
%           dVdv1=[dV/dv(1,derivIdx),dV/dv(2,derivIdx),dV/dv(3,derivIdx)].
%
%The volume of a tetrahedron can be expressed using a determinant as in
%[1]. The sum is implemented here using accurateSum. The sign of the
%determininat is retained if onlyPositive=false.
%
%EXAMPLE 1:
%As an example, we consider a regular tetrahedron, for which a simple exact
%expression for the volume is available as in [2]. The relative error
%bewteen this function and the exact solution is on the order of what one
%would expect due to finite precision limitations.
% x=sqrt(3)/3;
% h=sqrt(6)/3;
% d=sqrt(3)/6;
% v=[x,  -d,  -d, 0;
%    0, 1/2,-1/2, 0;
%    0,   0,   0, h];
% trueVol=1/(6*sqrt(2));
% vol=tetrahedronVolume(v,true);
% RelErr=(vol-trueVol)/trueVol
%
%EXAMPLE 2:
%In this example, we verify that the derivative is consistent with
%numerical differentiation. The absolute error is consistent with what one
%might expect due to finite precision limitations.
% v=randn(3,4);
% [~,dVdv1]=tetrahedronVolume(v);
% f=@(v1)tetrahedronVolume([v1,v(:,2:4)]);
% dVdv1NumDiff=numDiff(v(:,1),f,1);
% AbsErr=max(abs((dVdv1NumDiff-dVdv1)))
%
%REFERENCES:
%[1] Weisstein, Eric W. "Tetrahedron." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/Tetrahedron.html
%[2] Jackson, Frank and Weisstein, Eric W. "Regular Tetrahedron."
%    From MathWorld--A Wolfram Web Resource.
%    https://mathworld.wolfram.com/RegularTetrahedron.html 
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(derivIdx))
    derivIdx=1;
end

if(nargin<2||isempty(onlyPositive))
    onlyPositive=false;
end

x1=v(1,1);
x2=v(1,2);
x3=v(1,3);
x4=v(1,4);
y1=v(2,1);
y2=v(2,2);
y3=v(2,3);
y4=v(2,4);
z1=v(3,1);
z2=v(3,2);
z3=v(3,3);
z4=v(3,4);

V=accurateSum([-x3*y2*z1, x4*y2*z1, x2*y3*z1,-x4*y3*z1,-x2*y4*z1, x3*y4*z1,...
              x3*y1*z2,-x4*y1*z2,-x1*y3*z2, x4*y3*z2, x1*y4*z2,-x3*y4*z2,...
             -x2*y1*z3, x4*y1*z3, x1*y2*z3,-x4*y2*z3,-x1*y4*z3, x2*y4*z3,... 
              x2*y1*z4,-x3*y1*z4,-x1*y2*z4, x3*y2*z4, x1*y3*z4,-x2*y3*z4],true)/6;

%Make it positive.
if(onlyPositive)
    s=sign(V);
    V=s*V;
else
    s=1;
end
    
if(nargout>1)
    switch(derivIdx)
        case 1    
            dVdx=accurateSum([-y3*z2, y4*z2, y2*z3, -y4*z3, -y2*z4, y3*z4],true)/6;
            dVdy=accurateSum([x3*z2, -x4*z2, -x2*z3, x4*z3, x2*z4, -x3*z4],true)/6;
            dVdz=accurateSum([-x3*y2, x4*y2, x2*y3, -x4*y3, -x2*y4, x3*y4],true)/6;
        case 2
            dVdx=accurateSum([y3*z1, -y4*z1, -y1*z3, y4*z3, y1*z4, -y3*z4],true)/6;
            dVdy=accurateSum([-x3*z1, x4*z1, x1*z3, -x4*z3, -x1*z4, x3*z4],true)/6;
            dVdz=accurateSum([x3*y1, -x4*y1, -x1*y3, x4*y3, x1*y4, -x3*y4],true)/6;
        case 3
            dVdx=accurateSum([-y2*z1, y4*z1, y1*z2, -y4*z2, -y1*z4, y2*z4],true)/6;
            dVdy=accurateSum([x2*z1, -x4*z1, -x1*z2, x4*z2, x1*z4, -x2*z4],true)/6;
            dVdz=accurateSum([-x2*y1, x4*y1, x1*y2, -x4*y2, -x1*y4, x2*y4],true)/6;
        case 4
            dVdx=accurateSum([y2*z1, -y3*z1, -y1*z2, y3*z2, y1*z3, -y2*z3],true)/6;
            dVdy=accurateSum([-x2*z1, x3*z1, x1*z2, -x3*z2, -x1*z3, x2*z3],true)/6;
            dVdz=accurateSum([x2*y1, -x3*y1, -x1*y2, x3*y2, x1*y3, -x2*y3],true)/6;
        otherwise
            error('Invalid derivative index specified.')
    end
  
    dVdv=s*[dVdx,dVdy,dVdz];
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
