function [a,t,tLength,k]=BSplineInterpDerivMultiDim(t,tLength,a,k,numDerivs)
%%BSPLINEINTERPDERIVMULTIDIM Given a hypermatrix of multivariate b-spline
%            interpolation weights a and the collection of knot components
%            t, compute the modified weights and knots to interpolate the
%            multivariate derivative of the function being considered. the
%            function can be complex, but its arguments must be real.
%
%INPUTS: t The maxNumKnotsXnumDims set of knots fo r the interpolation
%          function. These are actually the values in each coordinate
%          dimension given in ascending order. The full knots are implied
%          [t1,t2,etc.]=ndgrid(t(1:tLengths(1),1),t(1:tLengths(2),2),...)
%          and the dimensions can be put together into the full set of
%          numDimsXtotalNumPoints points as tTotal=[t1(:).';t2(:).';...].
%          The first and last k(curDim)-1 knots in each dimension are
%          outside of the ends of the region with data or mark the ends of
%          the region with data. The number of rows is the maximum number
%          needed for all of the dimensions. tLength says how many items
%          in each column are actually used. These values must be real.
%  tLength A numDimX1 vector where tLength(i) says the number of elements
%          in t(:,i).
%        a A hypermatrix with numDim+1 indices containing the set of
%          coefficients for the b-splines covering all of the interpolation
%          intervals for all of the sets of interpolation values. If there
%          is only one set, the final dimension is unitary meaning that
%          there are effectively only numDim dimensions. The derivatives
%          will be applied to all sets present. These values can be real or
%          complex.
%        k A numDimsX1 or 1XnumDims set of the order of the b-splines in
%          each dimension. The value k-1 is the polynomial order of the
%          approximation being performed. If the order is the same in all
%          dimensions, then a single scalar can be passed.
% numDerivs A numDimsX1 or 1XnumDims vector where each dimension specified
%          the number of derivatives to take. Each element must be less
%          than k-1
%
%OUTPUTS: a The modified coefficient matrix. The dimensions are not
%           modified even though interpolation of the derivatives does not
%           require the extra elements.
%         t The modified knots.
%   tLength The modified vector indicating how many knows are present.
%         k The modified set of orders.
%
%One-dimensional derivatives are given by Equation 12b in Chapter X of [1].
%Chapter XVII describes the use of tensor product splines for more than one
%dimension. The generalization to more than one dimension comes from just
%evaluating the derivatives across each dimension one at a time.
%
%EXAMPLE 1:
%Given a function with 2D input that is fourth-order in x and y, we find
%fifth-order b-splines, assuring an exact fit. We then evaluate take the
%second derivative in x and then the third derivative in y. The
%interpolated and true results are plotted and can be seen to be the same.
% f=@(x,y)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y);
% %The second derivative of f with respect to x and the third with respect
% %to y.
% fDeriv2X3Y=@(x,y)24*(-4+12*x.^2).*y;
% numDims=2;
% numPointsX=10;
% numPointsY=11;
% tauLengths=[numPointsX;numPointsY];
% tau=zeros(max(tauLengths),numDims);
% tau(1:tauLengths(1),1)=linspace(-1.5,1.5,numPointsX);
% tau(1:tauLengths(2),2)=linspace(-1.5,1.5,numPointsY);
% %Note that using meshgrid would have put the elements in the wong order.
% %ndgrid must be used.
% [tau1,tau2]=ndgrid(tau(1:tauLengths(1),1),tau(1:tauLengths(2),2));
% y=f(tau1,tau2);
% 
% k=[5;5];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,y(:),k);
% numDerivs=[2;3];
% [a,t,tLength,k]=BSplineInterpDerivMultiDim(t,tLength,a,k,numDerivs);
% 
% numPointsX=100;
% numPointsY=101;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% [X,Y]=ndgrid(pointsX,pointsY);
% x=[X(:).';Y(:).'];
% 
% zDerivTrue=fDeriv2X3Y(X,Y);
% z=BSplineInterpValMultiDim(x,t,tLength,a,k);
% z=reshape(z,[numPointsX,numPointsY]);
% 
% figure(1)
% clf
% surface(X,Y,z,'EdgeColor','none')
% title('Interpolated Surface')
% figure(2)
% clf
% surface(X,Y,zDerivTrue,'EdgeColor','none')
% title('True Function')
% max(abs(zDerivTrue(:)-z(:)))
%One sees that the curves are the same and that the difference at any point
%is on the order of finite precision error. This is consistent with the
%fact that the order of the spline fit was chosen to be the same as the
%order of the true polynomial.
%
%EXAMPLE 2:
%This example is similar to example 1 except a 3D function is used. Here,
%we just compute the maximum error over a finer interpolation grid rather
%than plot anything. 
% f=@(x,y,z)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y).*(z.^4-3*z.^2+z);
% %The second derivative of f with respect to x, the third with respect to y,
% %and the first with respect to z.
% fDeriv2X3Y1Z=@(x,y,z)24*(-4+12*x.^2).*y.*(1-6*z+4*z.^3);
% numDims=3;
% numPointsX=10;
% numPointsY=11;
% numPointsZ=9;
% tauLengths=[numPointsX;numPointsY;numPointsZ];
% tau=zeros(max(tauLengths),numDims);
% tau(1:tauLengths(1),1)=linspace(-1.5,1.5,numPointsX);
% tau(1:tauLengths(2),2)=linspace(-1.5,1.5,numPointsY);
% tau(1:tauLengths(3),3)=linspace(-1.5,1.5,numPointsZ);
% %Note that using meshgrid would have put the elements in the wong order.
% %ndgrid must be used.
% [tau1,tau2,tau3]=ndgrid(tau(1:tauLengths(1),1),tau(1:tauLengths(2),2),tau(1:tauLengths(3),3));
% y=f(tau1,tau2,tau3);
% 
% k=[5;5;5];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,y(:),k);
% numDerivs=[2;3;1];
% [a,t,tLength,k]=BSplineInterpDerivMultiDim(t,tLength,a,k,numDerivs);
% 
% numPointsX=20;
% numPointsY=21;
% numPointsZ=19;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% pointsZ=linspace(-1.5,1.5,numPointsZ);
% [X,Y,Z]=ndgrid(pointsX,pointsY,pointsZ);
% x=[X(:).';Y(:).';Z(:).'];
% 
% zDerivTrue=fDeriv2X3Y1Z(X,Y,Z);
% z=BSplineInterpValMultiDim(x,t,tLength,a,k);
% max(abs(z(:)-zDerivTrue(:)))
%As was the case with Example 1, the maximum error is on the order of what
%one would expect due to finite precision limitations.
%
%EXAMPLE 3:
%This is like the first example, except multiple sets of points are used at
%once with multiple functions.
% f1=@(x,y)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y);
% f2=@(x,y)(x.^4-2*x.^3+x).*(y.^4-2*y.^2+y);
% %The second derivative of f1 and f2 with respect to x and the third with
% %respect to y.
% fDeriv2X3Y=@(x,y)[24*(-4+12*x.^2).*y,24*(-12*x+12*x.^2).*y];
% numDims=2;
% numPointsX=10;
% numPointsY=11;
% tauLengths=[numPointsX;numPointsY];
% tau=zeros(max(tauLengths),numDims);
% tau(1:tauLengths(1),1)=linspace(-1.5,1.5,numPointsX);
% tau(1:tauLengths(2),2)=linspace(-1.5,1.5,numPointsY);
% %Note that using meshgrid would have put the elements in the wong order.
% %ndgrid must be used.
% [tau1,tau2]=ndgrid(tau(1:tauLengths(1),1),tau(1:tauLengths(2),2));
% y=[f1(tau1(:),tau2(:)),f2(tau1(:),tau2(:))];
% 
% k=[5;5];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,y,k);
% numDerivs=[2;3];
% [a,t,tLength,k]=BSplineInterpDerivMultiDim(t,tLength,a,k,numDerivs);
% 
% numPointsX=100;
% numPointsY=101;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% [X,Y]=ndgrid(pointsX,pointsY);
% x=[X(:).';Y(:).'];
% 
% zDerivsTrue=reshape(fDeriv2X3Y(X(:),Y(:)),[numPointsX,numPointsY,2]);
% z=BSplineInterpValMultiDim(x,t,tLength,a,k);
% z=reshape(z,[numPointsX,numPointsY,2]);
% 
% figure(1)
% clf
% surface(X,Y,z(:,:,1),'EdgeColor','none')
% title('Interpolated Surface')
% figure(2)
% clf
% surface(X,Y,zDerivsTrue(:,:,1),'EdgeColor','none')
% title('True Function')
% max(max(abs(zDerivsTrue(:,:,1)-z(:,:,1))))
% 
% figure(3)
% clf
% surface(X,Y,z(:,:,2),'EdgeColor','none')
% title('Interpolated Surface')
% figure(4)
% clf
% surface(X,Y,zDerivsTrue(:,:,2),'EdgeColor','none')
% title('True Function')
% max(max(abs(zDerivsTrue(:,:,2)-z(:,:,2))))

%EXAMPLE 4:
%This example demonstrates differentiation of a bivariate complex function
%with respect to both of its real inputs.
% numPoints=100;
% xMin=0;
% xMax=3;
% yMin=0;
% yMax=3;
% x=linspace(xMin,xMax,numPoints);
% y=linspace(yMin,yMax,numPoints);
% [X,Y]=ndgrid(x,y);
% f=@(x,y)((sin(5*y)+1j*cos(10*x)).*exp(1j*(x.^2-2*y.^2)));
% d2fdxdy=@(x,y)(2*exp(1j*(x.^2-2*y.^2)).*(4*1j*x.*y.*cos(10*x)+5*1j*x.*cos(5*y)+4*y.*(-5*sin(10*x)+x.*sin(5*y))));
% 
% fxy=f(X(:),Y(:));
% k=[5;5];
% tau=[x(:),y(:)];
% tauLengths=[numPoints,numPoints];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,fxy,k);
% numDerivs=[1;1];
% [a,t,tLength,k]=BSplineInterpDerivMultiDim(t,tLength,a,k,numDerivs);
% 
% %Interpolate the derivative
% numPoints=100;
% x=linspace(xMin,xMax,numPoints);
% y=linspace(xMin,xMax,numPoints);
% [X,Y]=ndgrid(x,y);
% 
% fxy=d2fdxdy(X(:),Y(:));
% pts=[X(:).';Y(:).'];
% fxyInterp=reshape(BSplineInterpValMultiDim(pts,t,tLength,a,k),[numPoints,numPoints]);
% 
% fxy=reshape(fxy,[numPoints,numPoints]);
% 
% figure(1)
% clf
% hold on
% surface(X,Y,real(fxy),'EdgeColor','none')
% title('Real Derivative Values')
% colorbar()
% 
% figure(2)
% clf
% hold on
% surface(X,Y,imag(fxy),'EdgeColor','none')
% title('Imaginary Derivative Values')
% colorbar()
% 
% figure(3)
% clf
% hold on
% surface(X,Y,abs(real(fxyInterp)-real(fxy)),'EdgeColor','none')
% title('Real Interpolation Error')
% colorbar()
% 
% figure(4)
% clf
% hold on
% surface(X,Y,abs(imag(fxyInterp)-imag(fxy)),'EdgeColor','none')
% title('Imaginary Interpolation Error')
% colorbar()
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=length(tLength);
numA=size(a);

if(length(numA)==numDims)
    numSets=1;
else
    numSets=numA(end);
end

if(numDims==1)
    %Just use the 1D function if 1D inputs are passed.
    [t,a,k]=BSplineInterpDeriv(t(1:tLength),a,k,numDerivs);
    tLength=tLength-2*numDerivs;
    return; 
end

%When taking derivatives in each dimension, we have to go through all
%tuples of the other coordinates for all of the sets
for curDerivDim=1:numDims
    %If there are no derivatives in this dimension.
    if(numDerivs(curDerivDim)==0)
       continue;
    end
    
    %Next, we go through all tuples of values for the dimensions other than
    %this one. In each instance, we must take the derivative of the a terms
    %that are present. These indices select all of the dimensions that are
    %not the current one across which derivatives are being taken.
    for curSet=1:numSets
        idxList=[1:(curDerivDim-1),(curDerivDim+1):numDims];

        maxVals=numA(idxList)-1;
        numACur=numA(curDerivDim);
        curTuple=getNextTuple(numDims-1);
        idxCell=cell(1,numDims+1);
        idxCell{numDims+1}=curSet;
        tCur=t(1:tLength(curDerivDim),curDerivDim);
        kCur=k(curDerivDim);
        numDerivsCur=numDerivs(curDerivDim);

        outputIdx=1:(numACur-numDerivsCur);
        while(~isempty(curTuple))
            idxCell{curDerivDim}=':';%Mark the free dimension.

            %We have to select the elements in a to work on based on the
            %current tuple.
            for curIdx=1:(numDims-1)
                idxCell{idxList(curIdx)}=curTuple(curIdx)+1;
            end
            aCur=reshape(a(idxCell{:}),[numACur,1]);
            [~,aCur]=BSplineInterpDeriv(tCur,aCur,kCur,numDerivsCur);
            idxCell{curDerivDim}=outputIdx;
            a(idxCell{:})=aCur;

            curTuple=getNextTuple(curTuple,maxVals);
        end
    end
    
    %Update t and k for the current dimension.
    tLengthNew=tLength(curDerivDim)-2*numDerivsCur;
    t(1:tLengthNew,curDerivDim)=t((1+numDerivsCur):(tLength(curDerivDim)-numDerivsCur),curDerivDim);
    tLength(curDerivDim)=tLengthNew;
    k(curDerivDim)=k(curDerivDim)-numDerivsCur;
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
