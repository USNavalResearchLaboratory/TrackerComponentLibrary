function J=polarU2DCrossGrad(uList,systemType)
%%POLARU2DCROSSGRAD Given the direction cosine value u in 2D, obtain the
%              derivative of the polar azimuth angle with respect to u.
%              This is the derivative of the output of u2PolAng2D with
%              respect to its input. Note that this does not depend on the
%              rotation matrices MP and MU that are needed for the
%              u2PolAng2D function.
%
%INPUTS: uList A 1XnumPoints (for only u) or a 2XnumPoints (if full unit
%              vectors are given) set of direction cosines in 2D.
%   systemType An optional parameter specifying the axis from which the
%              angles are measured. Possible values are
%              0 (The default if omitted or an empty matrix is passed) The
%                azimuth angle is counterclockwise from the x axis.
%              1 The azimuth angle is measured clockwise from the y axis.
%
%OUTPUTS: J A 1XnumPoints or 2XnumPoints (if uList is 2XnumPoints) set of
%           derivatives of the azimuth angle with respect to u (first row)
%           and v (second row) evaluated at the given points.
%
%Note that u and v are not independent. Thus, in the comparison to finite
%differencing below (Example 2), when we offset u, we change v and vice
%versa.
%
%EXAMPLE 1:
%Here, we verify that the derivatives returned by this function are about
%equal to those returned via numeric differentiation (forward
%differencing).
% points=[0.1,0.2,-0.2,0,-0.9];%u points
% systemType=1;
% epsVal=1e-8;
% MP=randRotMat(2);
% MU=randRotMat(2);
% 
% az=u2PolAng2D(points,systemType,MP,MU);
% az1=u2PolAng2D(points+epsVal,systemType,MP,MU);
% JNumDiff=(az1-az)/epsVal;
% J=polarU2DCrossGrad(points,systemType);
% max(abs(JNumDiff-J))
%One will see that the difference is O(1e-8), which is a good agreement.
%
%EXAMPLE 2:
%This is the same as the first example, except both a u and a v component
%are provided, so derivatives with respect to each term are provided. The
%relative errors are O(1e-6) or O(1-e7), which is good agreement.
% points=[0.1,0.2,-0.2,-0.1,-0.9];%u
% points(2,:)=sqrt(1-points(1,:).^2);%v
% %Do the same thing, but switch the sign of v.
% points=[points,[points(1,:);-points(2,:)]];
% systemType=0;
% epsVal=1e-9;
% MP=randRotMat(2);
% MU=randRotMat(2);
% 
% az=u2PolAng2D(points,systemType,MP,MU);
% zDiff(1,:)=points(1,:)+epsVal;
% zDiff(2,:)=sign(points(2,:)).*sqrt(1-zDiff(1,:).^2);
% az1=u2PolAng2D(zDiff,systemType,MP,MU);
% JNumDiffU=(az1-az)/epsVal;
% zDiff(2,:)=points(2,:)+epsVal;
% zDiff(1,:)=sign(points(1,:)).*sqrt(1-zDiff(2,:).^2);
% az1=u2PolAng2D(zDiff,systemType,MP,MU);
% JNumDiffV=(az1-az)/epsVal;
% JNumDiff=[JNumDiffU;JNumDiffV];
% J=polarU2DCrossGrad(points,systemType);
% RelErr=max(max(abs(JNumDiff-J)))
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
    systemType=0; 
end

hasV=size(uList,1)>1;

N=size(uList,2);

J=zeros(1+hasV,N);
for curPoint=1:N
    u=uList(1,curPoint);
    
    if(hasV)
        v=uList(2,curPoint);
    else
        v=sqrt(1-u^2);
    end

    switch(systemType)
        case 0
            J(1,curPoint)=-1./v;
            if(hasV)
                J(2,curPoint)=1./u;
            end
        case 1
            J(1,curPoint)=1./v;
            if(hasV)
                J(2,curPoint)=-1./u;
            end
        otherwise
            error('Invalid system type specified.')
    end
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
