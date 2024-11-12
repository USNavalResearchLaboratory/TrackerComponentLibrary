function H=LorentzTransHessian(xTimePos,v,vectorType,c)
%%LORENTZTRANSHESSIAN This function evaluates the Hessian (matrix of
%     second partial derivatives) of the output of the LorentzTrans
%     function with respect to xTimePos and v.
%
%INPUTS: xTimePos The 4XN [time offset;position offset] values relative to
%            a common event, as viewed by the reference inertial coordinate
%            system. The exact format depends on the vectorType input.
%          v The 3X1 velocity vector of the origin of the second inertial
%            coordinate system measured with respect to the first inertial
%            coordinate system. Note that norm(vVec)<c, where c is the
%            speed of light in a vacuum.
% vectorType A string specifying the type of Lorentz transform matrix to
%            obtain. This can be
%            'Real' The default if this parameter is omitted. All of the
%                   entries in the Lorentz transforms matrix are real and
%                   the interval being transformed is assumed to have the
%                   real form xTimePos=[c*Delta t;Delta x']'
%            'RealAsymmetric' All of the entries in the Lorentz transforms
%                   matrix are real and the interval being transformed is
%                   assumed to have the real form
%                   xTimePos=[Delta t;Delta x']'
%          c The speed of light in a vacuum. If omitted or an empty matrix
%            is passed,t he default of c=Constants.speedOfLight, which has
%            united of meters per second, is used.
%
%OUTPUTS: H The 7X7X4 matrix of partial second derivatives, where H(:,:,i)
%           selects the ith dimensions of the output of LorentzTrans and
%           the rows and columns select the elements of the transformation.
%           Specifically, the elements are:
%           [d^2/(dtdt),   d^2/(dtdx),   d^2/(dtdy),   d^2/(dtdz),   d^2/(dtdxDot),   d^2/(dtdyDot),   d^2/(dtdzDot);
%            d^2/(dt),     d^2/(dxdx),   d^2/(dxdy),   d^2/(dxdz),   d^2/(dxdxDot),   d^2/(dxdyDot),   d^2/(dxdzDot);
%            d^2/(dydt),   d^2/(dydx),   d^2/(dydy),   d^2/(dydz),   d^2/(dydxDot),   d^2/(dydyDot),   d^2/(dydzDot);
%            d^2/(dzdt),   d^2/(dzdx),   d^2/(dzdy),   d^2/(dzdz),   d^2/(dzdxDot),   d^2/(dzdyDot),   d^2/(dzdzDot);
%            d^2/(dxDotdt),d^2/(dxDotdx),d^2/(dxDotdy),d^2/(dxDotdz),d^2/(dxDotdxDot),d^2/(dxDotdyDot),d^2/(dxDotdzDot);
%            d^2/(dyDotdt),d^2/(dyDotdx),d^2/(dyDotdy),d^2/(dyDotdz),d^2/(dyDotdxDot),d^2/(dyDotdyDot),d^2/(dyDotdzDot);
%            d^2/(dzDotdt),d^2/(dzDotdx),d^2/(dzDotdy),d^2/(dzDotdz),d^2/(dzDotdxDot),d^2/(dzDotdyDot),d^2/(dzDotdzDot)];
%
%This function implements explicit second derivatives of the function
%LorentzTrans. Expressions are given in [1].
%
%EXAMPLE:
%The results are compared to finite differencing. The results agree to more
%than 7 digits, so finite differencing is consistent with the explicit
%solution.
% c=Constants.speedOfLight;
% u=[1;-2;3];
% u=u/norm(u);
% v=0.9999*c*u;%A very fast velocity.
% vectorType='Real';
% xTimePos=[1e-3*c;1e3;-12e3;3e3];
% f=@(x)vec(LorentzTransGradient(x(1:4),x(5:7),vectorType,c));
% xHessNumDiff=permute(reshape(numDiff([xTimePos;v],f,28,6),[4,7,7]),[2,3,1]);
% H=LorentzTransHessian(xTimePos,v,vectorType,c);
% RelErr=max(abs((H(:)-xHessNumDiff(:))./xHessNumDiff(:)))
%
%REFERENCES:
%[1] D. F. Crouse, "Debiasing nonlinear transformations involving
%    correlated measurement components," in Proceedings of the 27th
%    International Conference on Infromation Fusion, Venice, Italy, 7-11
%    Jul. 2024.
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<3||isempty(vectorType))
    vectorType='Real';
end

LGrad=LorentzTransMatrixGrad(v,vectorType,c);
LHessian=LorentzTransMatrixHessian(v,vectorType,c);

H=zeros(7,7,4);
eVecI=zeros(4,1);
eVecJ=zeros(4,1);
for i=1:7
    if(i<5)
        eVecI(i)=1;
    end

    for j=i:7
        if(i<5&&j<5)
            %Second derivative involving only position compoennts are zero.
            continue;
        end

        if(j<5)
            eVecJ(j)=1;
            H(i,j,:)=reshape(LGrad(:,:,i-4)*eVecJ,[1,1,4]);
            H(j,i,:)=H(i,j,:);
            eVecJ(j)=0;
        elseif(i<5)
            H(i,j,:)=reshape(LGrad(:,:,j-4)*eVecI,[1,1,4]);
            H(j,i,:)=H(i,j,:);
        else%Both i and j are >5 and thus represent velocities.
            H(i,j,:)=reshape(LHessian(:,:,i-4,j-4)*xTimePos,[1,1,4]);
            H(j,i,:)=H(i,j,:);
        end
    end

    if(i<5)
        eVecI(i)=0;
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
