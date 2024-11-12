function xvGrad=LorentzTransGradient(xTimePos,v,vectorType,c)
%%LORENTZTRANSGRADIENT This function evaluates the gradient (matrix of
%     first partial derivatives) of the output of the LorentzTrans function
%     with respect to xTimePos and v.
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
%            'Real' (The default if this parameter is omitted or an empty
%                   matrix is passed). All of the entries in the Lorentz
%                   transform matrix are real and the interval being
%                   transformed is assumed to have the real form
%                   xTimePos=[c*Delta t;Delta x']'
%            'RealAsymmetric' All of the entries in the Lorentz transform
%                   matrix are real and the interval being transformed is
%                   assumed to have the real form
%                   xTimePos=[Delta t;Delta x']'
%          c The speed of light in a vacuum. If omitted or an empty matrix
%            is passed,t he default of c=Constants.speedOfLight, which has
%            united of meters per second, is used.
%
%OUTPUTS: xvGrad A 4X7 matrix where the rows correspond to the elements of
%                the Lorentz transformed xTimePos and the columns indicate
%                what derivative is taken, with the elements of xTimePos
%                being the first 4 columns and the element of v being the
%                final 3 columns.
%
%This function implements explicit first derivatives of the function
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
% vectorType='RealAsymmetric';
% xTimePos=[1e-3;1e3;-12e3;3e3];
% f=@(x)LorentzTrans(x(1:4),x(5:7),vectorType,c);
% xGradNumDiff=reshape(numDiff([xTimePos;v],f,4,4),[4,7]);
% xGrad=LorentzTransGradient(xTimePos,v,vectorType,c);
% RelErr=max(abs((xGrad(:)-xGradNumDiff(:))./xGradNumDiff(:)))
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

L=LorentzTransMatrix(v,vectorType,c);
LGrad=LorentzTransMatrixGrad(v,vectorType,c);

xvGrad=zeros(4,7);

%The non-velocity gradient contribution.
xvGrad(1:4,1:4)=L;
%The contribution due to the velocity components.
xvGrad(1:4,5)=LGrad(:,:,1)*xTimePos;
xvGrad(1:4,6)=LGrad(:,:,2)*xTimePos;
xvGrad(1:4,7)=LGrad(:,:,3)*xTimePos;

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
