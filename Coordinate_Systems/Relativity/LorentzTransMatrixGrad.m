function LGrad=LorentzTransMatrixGrad(v,vectorType,c)
%%LORENTZTRANSMATRIXGRAD Evaluate the first derivatives of the special
%   relativity Lorentz transformation matrix from the function
%   LorentzTransMatrix with respect to the velocity vector of the second
%   inertial coordinate system.
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
%            is passed, the default of c=Constants.speedOfLight, which has
%            united of meters per second, is used.
%
%OUTPUT: LGrad A 4X4X3 matrix such that LGrad(:,:,k) is the gradient of the
%              Lorentz transformation matrix with respect to v(k).
%
%This function just implements the analytic first derivatives of the real
%transformation matrices in LorentzTransMat.
%
%EXAMPLE:
%The results are compared to finite differencing. The results agree to more
%than 9 digits, so finite differencing is consistent with the explicit
%solution.
% c=Constants.speedOfLight;%Assumed propagation speed.
% u=[1;-2;3];
% u=u/norm(u);
% v=0.999*c*u;%A very fast velocity.
% vectorType='RealAsymmetric';
% f=@(x)vec(LorentzTransMatrix(x,vectorType));
% LGradNumDiff=reshape(numDiff(v,f,16,3),[4,4,3]);
% LGrad=LorentzTransMatrixGrad(v,vectorType);
% RelErr=max(abs((LGrad(:)-LGradNumDiff(:))./LGradNumDiff(:)))
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
v2=v'*v;
v4=v2*v2;

gamma=1/sqrt(1-v2/c2);
gradGamma=c*v/(c2-v2)^(3/2);

LGrad=zeros(4,4,3);
ek=zeros(3,1);
switch(vectorType)
    case 'Real'
        for k=1:3
            ek(k)=1;
        
            LGrad(1,1,k)=gradGamma(k);
            LGrad(1,2:4,k)=-gradGamma(k)/c*v'-gamma/c*ek';
            LGrad(2:4,1,k)=LGrad(1,2:4,k)';
            LGrad(2:4,2:4,k)=gradGamma(k)*vOuter/v2+(gamma-1)*((ek*v'+v*ek')/v2-2*v(k)*vOuter/v4);
            ek(k)=0;
        end
    case 'RealAsymmetric'
        for k=1:3
            ek(k)=1;
        
            LGrad(1,1,k)=gradGamma(k);
            LGrad(1,2:4,k)=-gradGamma(k)/c2*v'-gamma/c2*ek';
            LGrad(2:4,1,k)=-gradGamma(k)*v-gamma*ek;
            LGrad(2:4,2:4,k)=gradGamma(k)*vOuter/v2+(gamma-1)*((ek*v'+v*ek')/v2-2*v(k)*vOuter/v4);
            ek(k)=0;
        end
    otherwise
        error('Invalid vectorType given');
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
