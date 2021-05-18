function [a,t,tLength,k]=BSplineInterpIntMultiDim(t,tLength,a,k,intDims)
%%BSPLINEINTERPINTMULTIDIM Given a hypermatrix of multivariate b-spline
%            interpolation weights a and a set of knots t,compute the
%            modified weights and knots to interpolate the integral of the
%            a given point using b-spline interpolation. If all knots at
%            the ith boundary are repeated k(i) times, then an integral
%            over the ith dimension will start from the beginning of the
%            interpolation region. Otherwise, this can be considered an
%            indefinite integral with a particular additive constant and
%            differencing and be used to get a definite integral.
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
%  intDims This is a numDimsX1 or a 1XnumDims boolean vector indicating
%          over which dimensions the integrals should be taken.
%
%OUTPUTS: a The modified coefficient matrix. This is one larger in each
%           dimensions where an integral was taken.
%         t The modified knots.
%   tLength The modified vector indicating how many knows are present.
%         k The modified set of orders.
%
%Equation 22 in Chapter X of [1] has 1D integration. Chapter XVII describes
%the use of tensor product splines for more than one dimension. The
%generalization to more than one dimension comes from just evaluating the
%integrals across each dimension one at a time, if they are taken.
%
%EXAMPLE 1:
% f=@(x,y,z)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y).*(z.^4-3*z.^2+z);
% %The definite integral of f with respect to x and z starting at -1.5 is
% fIntXZ=@(x,y,z)((-891+240*x.^2-320*x.^3+96*x.^5).*y.*(1-2*y+y.^3).*(-477+80*z.^2-160*z.^3+32*z.^5))/76800;
% 
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
% intDims=[1;0;1];
% [a,t,tLength,k]=BSplineInterpIntMultiDim(t,tLength,a,k,intDims);
% 
% numPointsY=21;
% numPointsZ=19;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% pointsZ=linspace(-1.5,1.5,numPointsZ);
% [X,Y,Z]=ndgrid(pointsX,pointsY,pointsZ);
% x=[X(:).';Y(:).';Z(:).'];
% 
% zIntTrue=fIntXZ(X,Y,Z);
% z=BSplineInterpValMultiDim(x,t,tLength,a,k);
% max(abs(z(:)-zIntTrue(:)))
%
%EXAMPLE 2:
%This is similar to example 1, except two sets of values are fitted,
%integrated, and interpolated at once.
% f1=@(x,y,z)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y).*(z.^4-3*z.^2+z);
% f2=@(x,y,z)(x.^3-2*x.^2+x).*(y.^4-2*y.^2+y).*(2*z.^3-2*z.^2+1);
% %The definite integral of f with respect to x and z starting at -1.5 is
% fIntXZ1=@(x,y,z)((-891+240*x.^2-320*x.^3+96*x.^5).*y.*(1-2*y+y.^3).*(-477+80*z.^2-160*z.^3+32*z.^5))/76800;
% fIntXZ2=@(x,y,z)(((-891+16*x.^2.*(6+x.*(-8+3*x))).*y.*(1-2*y+y.^3).*(-315+16*z.*(6+z.^2.*(-4+3*z))))/18432);
% 
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
% y1=f1(tau1,tau2,tau3);
% y2=f2(tau1,tau2,tau3);
% y=[y1(:),y2(:)];
% 
% k=[5;5;5];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,y,k);
% intDims=[1;0;1];
% [a,t,tLength,k]=BSplineInterpIntMultiDim(t,tLength,a,k,intDims);
% 
% numPointsY=21;
% numPointsZ=19;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% pointsZ=linspace(-1.5,1.5,numPointsZ);
% [X,Y,Z]=ndgrid(pointsX,pointsY,pointsZ);
% x=[X(:).';Y(:).';Z(:).'];
% 
% zIntTrue1=fIntXZ1(X,Y,Z);
% zIntTrue2=fIntXZ2(X,Y,Z);
% z=BSplineInterpValMultiDim(x,t,tLength,a,k);
% max(abs(z(:,1)-zIntTrue1(:)))
% max(abs(z(:,2)-zIntTrue2(:)))
%
%EXAMPLE 3:
%This is an example of a bivariate complex function whose integral over the
%first dimension is interpolated.numPoints=80;
% xMin=0;
% xMax=3;
% yMin=0;
% yMax=3;
% x=linspace(xMin,xMax,numPoints);
% y=linspace(yMin,yMax,numPoints);
% [X,Y]=ndgrid(x,y);
% f=@(x,y)(sin(5*y)+1j*cos(10*x));
% %The integral from 0 in the x dimension.
% fIntx=@(x,y)((1/10)*1j*sin(10*x)+x.*sin(5*y));
% 
% fxy=f(X(:),Y(:));
% k=[5;5];
% tau=[x(:),y(:)];
% tauLengths=[numPoints,numPoints];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,fxy,k);
% intDims=[1;0];
% [a,t,tLength,k]=BSplineInterpIntMultiDim(t,tLength,a,k,intDims);
% 
% %Interpolate the derivative
% numPoints=100;
% x=linspace(xMin,xMax,numPoints);
% y=linspace(xMin,xMax,numPoints);
% [X,Y]=ndgrid(x,y);
% 
% fxy=fIntx(X(:),Y(:));
% pts=[X(:).';Y(:).'];
% fxyInterp=reshape(BSplineInterpValMultiDim(pts,t,tLength,a,k),[numPoints,numPoints]);
% fxy=reshape(fxy,[numPoints,numPoints]);
% 
% figure(1)
% clf
% hold on
% surface(X,Y,real(fxy),'EdgeColor','none')
% title('Real Integral Values')
% colorbar()
% 
% figure(2)
% clf
% hold on
% surface(X,Y,imag(fxy),'EdgeColor','none')
% title('Imaginary Integral Values')
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
    numA=[numA,1];
else
    numSets=numA(end);
end

tNew=zeros(max(tLength)+2,numDims);
for curDim=1:numDims
    if(intDims(curDim)==0)
        tNew(1:tLength(curDim),curDim)=t(1:tLength(curDim),curDim);
        continue;
    end

    temp=zeros(1,numDims);
    temp(curDim)=1;
    aNew=zeros(numA+[temp,0]);
    
    %Next, we go through all tuples of values for the dimensions other than
    %this one. In each instance, we must take the integral of the a terms
    %that are present. These indices select all of the dimensions that are
    %not the current one across which derivatives are being taken.
    for curSet=1:numSets
        idxList=[1:(curDim-1),(curDim+1):numDims];

        maxVals=numA(idxList)-1;
        numACur=numA(curDim);
        curTuple=getNextTuple(numDims-1);
        idxCell=cell(1,numDims+1);
        idxCell{numDims+1}=curSet;
        tCur=t(1:tLength(curDim),curDim);
        kCur=k(curDim);

        outputIdx=1:(numACur+1);
        while(~isempty(curTuple))
            idxCell{curDim}=1:numA(curDim);%Mark the free dimension.

            %We have to select the elements in a to work on based on the
            %current tuple.
            for curIdx=1:(numDims-1)
                idxCell{idxList(curIdx)}=curTuple(curIdx)+1;
            end

            aCur=reshape(a(idxCell{:}),[numACur,1]);

            [~,aCur]=BSplineInterpInt(tCur,aCur,kCur);
            idxCell{curDim}=outputIdx;
            aNew(idxCell{:})=aCur;

            curTuple=getNextTuple(curTuple,maxVals);
        end
    end
    numA(curDim)=numA(curDim)+1;
    a=aNew;
    k(curDim)=k(curDim)+1;
    tNew(1:(tLength(curDim)+2),curDim)=[t(1,curDim);t(1:tLength(curDim),curDim);t(tLength(curDim),curDim)];
    tLength(curDim)=tLength(curDim)+2;
end
a=aNew;
t=tNew;

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
