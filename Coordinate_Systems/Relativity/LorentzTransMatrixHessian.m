function LHessian=LorentzTransMatrixHessian(v,vectorType,c)
%%LORENTZTRANSMATRIXHESSIAN Evaluate the second derivatives of the special
%   relativity Lorentz transformation matrix from the function
%   LorentzTransMatrix with respect to the elements of the velocity vector
%   of the second inertial coordinate system.
%
%INPUTS:v The 3X1 velocity vector of the origin of the second inertial
%         coordinate system measured with respect to the first inertial
%         coordinate system. Note that norm(vVec)<=c, where c is the
%         speed of light.
% vectorType A string specifying the type of Lorentz transform matrix to
%            obtain. This can be
%            'Real' (The default if this parameter is omitted or an empty
%                   matrix is passed). All of the entries in the Lorentz
%                   transform matrix are real and the interval being
%                   transformed is assumed to have the real form
%                   z=[c*Delta t;Delta x']'
%            'RealAsymmetric' All of the entries in the Lorentz transform
%                   matrix are real and the interval being transformed is
%                   assumed to have the real form z=[Delta t;Delta x']'
%          c The speed of light in a vacuum. If omitted or an empty matrix
%            is passed,t he default of c=Constants.speedOfLight, which has
%            united of meters per second, is used.
%
%OUTPUTS: LHessian A 4X4X3X3 matrix such that LHEssian(:,:,i,j) is the
%                  second derivative of the Lorentz transform matrix with
%                  respect to v(i) and v(j).
%
%This function just implements the analytic second derivatives of the real
%transformation matrices in LorentzTransMat.
%
%EXAMPLE:
%The results are compared to finite differencing. The results agree to more
%than 10 digits, so finite differencing is consistent with the explicit
%solution.
% c=Constants.speedOfLight;%Assumed propagation speed.
% u=[1;-2;3];
% u=u/norm(u);
% v=0.999*c*u;%A very fast velocity.
% vectorType='RealAsymmetric';
% 
% H=LorentzTransMatrixHessian(v,vectorType,c);
% f=@(x)vec(LorentzTransMatrixGrad(x,vectorType,c));
% LHessianNumDiff=reshape(numDiff(v,f,48,3),[4,4,3,3]);
% RelErr=max(abs((H(:)-LHessianNumDiff(:))./LHessianNumDiff(:)))
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<2||isempty(vectorType))
    vectorType='Real';
end

c2=c*c;

vOuter=v*v';
vEls2=v.*v;
v2=sum(vEls2);
v4=v2*v2;
v6=v4*v2;

gamma=1/sqrt(1-v2/c2);
gradGamma=c*v/((c2-v2)^(3/2));

gravVTerms=zeros(3,3,3);
ei=zeros(3,1);
for i=1:3
    ei(i)=1;
    %This is the first derivative matrix of (v*v')/norm(v)^2 with respect
    %to v(i).
    gravVTerms(:,:,i)=(ei*v'+v*ei')/v2-2*v(i)*vOuter/v4;
    ei(i)=0;
end

denomTerm=(c^2-v2)^(5/2);
gammaHessian=zeros(3,3);
for k1=1:3
    gammaHessian(k1,k1)=c*(c2+3*vEls2(k1)-v2)/(denomTerm);

    for k2=(k1+1):3
        gammaHessian(k1,k2)=3*c*v(k1)*v(k2)/denomTerm;
        gammaHessian(k2,k1)=gammaHessian(k1,k2);
    end
end

LHessian=zeros(4,4,3);
ei=zeros(3,1);
ej=zeros(3,1);
for i=1:3
    ei(i)=1;

    LHessian(1,1,i,i)=gammaHessian(i,i);
    LHessian(1,2:4,i,i)=-gammaHessian(i,i)/c2*v'-2*gradGamma(i)/c^2*ei';
    LHessian(2:4,1,i,i)=-gammaHessian(i,i)*v-2*gradGamma(i)*ei;

    %This is the second derivative matrix of v*v'/norm(v)^2 with respect to
    %v(i).
    derivTerm2=2*(ei*ei')/v2-4*v(i)*(ei*v'+v*ei')/v4-2*(vOuter)/v4+8*v(i)^2*(vOuter)/v6;
    LHessian(2:4,2:4,i,i)=gammaHessian(i,i)*vOuter/v2+2*gradGamma(i)*gravVTerms(:,:,i)+(gamma-1)*derivTerm2;

    for j=(i+1):3
        ej(j)=1;

        %This is the second derivative of v*v'/norm(v)^2 with respect to
        %v(i) and v(j) with i~=j.
        derivTerm2=(ei*ej'+ej*ei')/v2-2*v(j)*(ei*v'+v*ei')/v4-2*v(i)*(ej*v'+v*ej')/v4+8*v(i)*v(j)*(v*v')/v6;

        LHessian(1,1,i,j)=gammaHessian(i,j);
        LHessian(1,2:4,i,j)=-gammaHessian(i,j)/c2*v'-gradGamma(i)/c^2*ej'-gradGamma(j)/c^2*ei';
        LHessian(2:4,1,i,j)=-gammaHessian(i,j)*v-gradGamma(i)*ej-gradGamma(j)*ei;
        LHessian(2:4,2:4,i,j)=gammaHessian(i,j)*vOuter/v2+gradGamma(i)*gravVTerms(:,:,j)+gradGamma(j)*gravVTerms(:,:,i)+(gamma-1)*derivTerm2;
    
        LHessian(:,:,j,i)=LHessian(:,:,i,j);

        ej(j)=0;
    end
    ei(i)=0;
end

switch(vectorType)
    case 'Real'
        LHessian(1,2:4,:,:)=LHessian(1,2:4,:,:)*c;
        LHessian(2:4,1,:,:)=LHessian(2:4,1,:,:)/c;
    case 'RealAsymmetric'
    otherwise
        error('Invalid vectorType given');
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
