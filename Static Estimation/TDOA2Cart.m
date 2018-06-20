function [zCart,exitCode]=TDOA2Cart(TDOA,lRx1,lRx2,c,AbsTol,RelTol,params)
%%TDOA2CART Given the minimum number of time difference of arrival (TDOA)
%           measurements between sensors necessary to determine the
%           location of an object (2 in 2D and 3 in 3D), solve for the
%           location of the object by solving a set of simultaneous
%           multivariate polynomial equations. Thresholds can be specified
%           for discarding complex solutions. Unlike the function
%           TDOAOnlyStaticLocEst, this function can work with the minimum
%           number of values needed to mathematically solve the problem,
%           but this function cannot make use of additional measurements.
%
%INPUTS: TDOA A vector containing 2 (for 2D) or 3 (for 3D) TDOA
%             measurements. The TDOA measurements correspond to differences
%             taken between sensors whose locations are given in lRx2 and
%             lRx1.
%        lRx1 A 2X2 or 3X3 matrix of 2 (or 2D) or 3 (for 3D) positions of
%             the receiver locations that form the reference time from
%             which the TDOA measurement is computed. lRx1(:,i) is the
%             location of the ith receiver. If all of the reference times
%             are the same, then this can be a single 2X1 or 3X1 vector.
%        lRx2 A 2X2 or 3X3 matrix of 2 (or 2D) or 3 (for 3D) positions of
%             the receiver locations from whose reception time the time of
%             arrival at lRx1 is subtracted to get the TDOA measurements.
%           c The propagation speed in the medium in question. If this
%             parameter is omitted or an empty matrix is passed, the
%             default value of Constants.speedOfLight is used.
% AbsTol,RelTol In some geometries of targets and
%             recievers and when noise is added to the measurements,
%             meaningless complex solutions can be produced. AbsTol and
%             RelTol are tolerances for determining whether a solution is
%             numerically real and should be kept or numerically complex
%             and should be discarded. A solution x is kept if
%             abs(imag(x))<AbsTol or if abs(imag(x))<RelTol*abs(real(x)).
%             If omitted or empty matrices are passed, the default values
%             of AbsTol=1e-9 and RelTol=1e-7 are used. Additionally, AbsTol
%             and RelTol are used to eliminate false solutions that arise
%              due to the TDOAs being squared when
%             formulating the polynomial problem. These solutions are
%             eliminated if ~(absDiff<AbsTol ||
%             absDiff<RelTol*abs(TDOA(curEq))) is true. This criterion is
%             used rather than just comparing the signs of TDOA
%             computations to improve performance when TDOA values are
%             small/ zero.
%      params An optional structure of parameters that are for the
%             polyRootsMultiDim function, which is used by this function.
%             Possible members are maxDegIncreases and useMotzkinNull, both
%             of which are described in the function polyRootsMultiDim.
%             Default values for the function polyRootsMultiDim will be
%             used if this parameter is omitted or an empty matrix is
%             passed.
%
%OUTPUTS: zCart The 2XnumSol or 3XnumSol set of 2D or 3D real solutions for
%               the location of the target.
%      exitCode This is the exit code returned by the function
%               polyRootsMultiDim, which is used to solve the problem.
%
%The algorithm is based on the simultaneous multivariate polynomial
%formulation of [1]. Note that in 3D, the third measurement must not be a
%permutation of the other two, meaning at least four sensors are required.
%For example, one cannot use measurements of sensor1-sensors2,
%sensor1-sensor3 and sensor2-sensor3, since the last difference depends on
%the previous two. However, one could use sensor2-sensor4 and it would be
%independent.
%
%EXAMPLE 1:
%The example here is similar to the 2D example from [1]. Note that the
%TDOA numbers given in the example are incorrect. The correct numbers are
%used here.
% S1=[9;39];
% S2=[65;10];
% S3=[64;71];
% c=341;%The propagation speed.
% t=[27;42];%The true target location.
% 
% lRx1=zeros(2,2);
% lRx2=zeros(2,2);
% TDOA=zeros(2,1);
% %First TDOA pair is S2-S1
% lRx1(:,1)=S1;
% lRx2(:,1)=S2;
% TDOA(1)=(norm(t-S2)-norm(t-S1))/c;
% %Second TDOA pair is S3-S1
% lRx1(:,2)=S1;
% lRx2(:,2)=S3;
% TDOA(2)=(norm(t-S3)-norm(t-S1))/c;
% zCart=TDOA2Cart(TDOA,lRx1,lRx2,c)
%One will get the correct solution [27;42], within finite precision limits.
%
%EXAMPLE 2:
%Here, we give a 3D example.
% S1=[9;39;100];
% S2=[65;10;-60];
% S3=[64;71;43];
% S4=[-128;6;12];
% c=341;%The propagation speed.
% t=[27;0;-42];%The true target location.
% 
% lRx1=zeros(3,3);
% lRx2=zeros(3,3);
% TDOA=zeros(3,1);
% %First TDOA pair is S2-S1
% lRx1(:,1)=S1;
% lRx2(:,1)=S2;
% TDOA(1)=(norm(t-S2)-norm(t-S1))/c;
% %Second TDOA pair is S3-S1
% lRx1(:,2)=S1;
% lRx2(:,2)=S3;
% TDOA(2)=(norm(t-S3)-norm(t-S1))/c;
% %Third TDOA pair is S3-S4
% lRx1(:,3)=S4;
% lRx2(:,3)=S3;
% TDOA(3)=(norm(t-S3)-norm(t-S4))/c;
% zCart=TDOA2Cart(TDOA,lRx1,lRx2,c)
%One will get the correct solution [27;0;-42].
%
%REFERENCES:
%[1] M. P. Williams, "Solving polynomial equations using linear algebra,"
%    Johns Hopkins Technical Digest, vol. 28, no. 4, pp. 354-363, 2010.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<5||isempty(AbsTol))
    AbsTol=1e-9;
end

if(nargin<6||isempty(RelTol))
    RelTol=1e-7; 
end

if(nargin<7)
    params=[]; 
end

maxDegIncreases=[];
useMotzkinNull=[];
if(~isempty(params))
    if(isfield(params,'maxDegIncreases'))
        maxDegIncreases=params.maxDegIncreases; 
    end

    if(isfield(params,'useMotzkinNull'))
        useMotzkinNull=params.useMotzkinNull;
    end
end

numDim=size(lRx1,1);

%If all of the reference points are the same.
if(size(lRx1,2)==1)
    lRx1=repmat(lRx1,[1,numDim]);
end

%For the definition of lRx1 and lRx2 used in the equations, it is the
%opposite of the definition above, so we swap them.
temp=lRx1;
lRx1=lRx2;
lRx2=temp;

u=(lRx1+lRx2)/2;
v=(lRx1-lRx2)/2;

switch(numDim)
    case 2%The 2D Case
        xPolys=cell(2,1);
        for curMeas=1:2
            %The equation from the application section of [1] substituting
            %y=x-u (the sign of u was a typo in the paper), so we can
            %rewrite it in terms of x, is
            %-2*(u1*v1^2+u2*v1*v2-u1*delta^2)*x1+(v1^2-delta^2)*x1^2-2*(u1*v1*v2+u2*v2^2-u2*delta^2)*x2+2*v1*v2*x1*x2+(v2^2-delta^2)*x2^2+(-u1^2*delta^2-u2^2*delta^2-v1^2*delta^2-v2^2*delta^2+u1^2*v1^2+2*u1*u2*v1*v2+u2^2*v2^2+delta^4)==0            
            u1=u(1,curMeas);
            u2=u(2,curMeas);

            v1=v(1,curMeas);
            v2=v(2,curMeas);
            delta=c*TDOA(curMeas)/2;
            
            xPoly=zeros(3,3);
            xPoly(1+1,0+1)=-2*(u1*v1^2+u2*v1*v2-u1*delta^2);
            xPoly(2+1,0+1)=v1^2-delta^2;
            xPoly(0+1,1+1)=-2*(u1*v1*v2+u2*v2^2-u2*delta^2);
            xPoly(1+1,1+1)=2*v1*v2;
            xPoly(0+1,2+1)=v2^2-delta^2;
            xPoly(0+1,0+1)=-u1^2*delta^2-u2^2*delta^2-v1^2*delta^2-v2^2*delta^2+u1^2*v1^2+2*u1*u2*v1*v2+u2^2*v2^2+delta^4;
            xPolys{curMeas}=xPoly;
        end
    case 3%The 3D Case
        xPolys=cell(3,1);
        
        for curMeas=1:3
            %The equation from the application section of [1] substituting
            %y=x-u (the sign of u was a typo in the paper), so we can
            %rewrite it in terms of x, is
            %(u1*v1+u2*v2+u3*v3)^2+2*v1*v2*x1*x2+2*v1*v3*x1*x3+2*v2*v3*x2*x3-(u1^2+u2^2+u3^2+v1^2+v2^2+v3^2)*delta^2+delta^4+x1^2*(v1-delta)*(v1+delta)+x2^2*(v2-delta)*(v2+delta)+x3^2*(v3-delta)*(v3+delta)+x1*(-2*v1*(u1*v1+u2*v2+u3*v3)+2*u1*delta^2)+x2*(-2*v2*(u1*v1+u2*v2+u3*v3)+2*u2*delta^2)+x3*(-2*v3*(u1*v1+u2*v2+u3*v3)+2*u3*delta^2)
            u1=u(1,curMeas);
            u2=u(2,curMeas);
            u3=u(3,curMeas);

            v1=v(1,curMeas);
            v2=v(2,curMeas);
            v3=v(3,curMeas);
            delta=c*TDOA(curMeas)/2;
            
            xPoly=zeros(3,3,3);
            xPoly(1+1,0+1,0+1)=-2*v1*(u1*v1+u2*v2+u3*v3)+2*u1*delta^2;
            xPoly(2+1,0+1,0+1)=(v1-delta)*(v1+delta);
            xPoly(0+1,1+1,0+1)=-2*v2*(u1*v1+u2*v2+u3*v3)+2*u2*delta^2;
            xPoly(0+1,2+1,0+1)=(v2-delta)*(v2+delta);
            xPoly(1+1,1+1,0+1)=2*v1*v2;
            xPoly(0+1,0+1,1+1)=-2*v3*(u1*v1+u2*v2+u3*v3)+2*u3*delta^2;
            xPoly(1+1,0+1,1+1)=2*v1*v3;
            xPoly(0+1,1+1,1+1)=2*v2*v3;
            xPoly(0+1,0+1,2+1)=(v3-delta)*(v3+delta);
            xPoly(0+1,0+1,0+1)=(u1*v1+u2*v2+u3*v3)^2-(u1^2+u2^2+u3^2+v1^2+v2^2+v3^2)*delta^2+delta^4;
            xPolys{curMeas}=xPoly;
        end
    otherwise
        error('The dimensionality of the locations is invalid.')
end

%Plug the result into a multivariate polynomial solver.
[theRoots,exitCode]=polyRootsMultiDim(xPolys,maxDegIncreases,useMotzkinNull);

%Eliminate complex solutions.
sel=sum((abs(imag(theRoots))<AbsTol) | abs(imag(theRoots))<RelTol*abs(real(theRoots)),1)~=0;
zCart=real(theRoots(:,sel));

%Now, due to squaring in the equations, the sign of the solutions will
%generally be incorrect. Thus, we shall discard any solutions where the
%sign of the computed TDOA measurements does not match the sign of the
%actual TDOA measurements. However, this poses problems for TDOA values
%near zero. Thus, as opposed to comparing the sign, we use the same AbsTol
%and RelTol as for determining whether something is complex to determine
%whether we should discard a solution.
numSol=size(zCart,2);
sel=true(numSol,1);
for curSol=1:numSol
   for curEq=1:numDim
        TDOAComp=(norm(zCart(:,curSol)-lRx1(:,curEq))-norm(zCart(:,curSol)-lRx2(:,curEq)))/c;
        
        absDiff=abs(TDOAComp-TDOA(curEq));

        if(~(absDiff<AbsTol || absDiff<RelTol*abs(TDOA(curEq))))
            sel(curSol)=false;
            break;
        end
   end
end

%Make sure that we do not end up with no solutions.
if(all(sel==false))
    sel(1)=true;
end

zCart=zCart(:,sel);

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
