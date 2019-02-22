function [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,y,k,numCentralKnots,useSparse)
%%BSPLINEPOLYFITMULTIDIM Given a set of real vector (or scalar) points as
%         well as real or complex scalar function values, this function
%         returns the knots and coefficients for a set of b-splines that
%         can be used to interpolate to a desired polynomial degree within
%         the sampled region of the function. B-splines can have continuous
%         derivatives across the entire sampled region and can be easily
%         differentiated. This function determines the knots to use by
%         taking the tensor product of the approximate knots from
%         Equation 10 of Chapter 13 of [1]. The knots are not optimal. The
%         first and last knots in each dimensions are repeated. This
%         function can provide interpolation coefficients for multiple sets
%         of values at once. 
%
%INPUTS: tau A numPointsXnumDims matrix of the coordinates of the points in
%            each dimension in ascending order. The full grid of points
%            in all dimensions could be found using
%            [tau1,tau2,etc.]=ndgrid(tau(1:tauLengths(1),1),tau(1:tauLengths(2),2),...)
%            and the dimensions can be put together into the full set of
%            numDimsXtotalNumPoints points as
%            tauTotal=[tau1(:).';tau2(:).';...]. Note that using meshgrid
%            in 2D/3D will put the points in the wrong order; one has to
%            use the ordering of ndgrid. These values must be real.
% tauLengths A numDimsX1 or a 1XnumDims vector this lists the number of
%            points in each dimension of tau. The total number of sample
%            points is totalNumPoints=prod(tauLengths).
%          y This is the totalNumPointsXnumSets full sets of function
%            values taken at all of the points implied by the vectors in
%            tau for numSets different interpolation problems. The ordering
%            of the points is the same as the ordering of tauTotal, which
%            is described above. These values can be real or complex.
%          k A numDimsX1 or 1XnumDims set of the order of the b-splines in
%            each dimension. The value k-1 is the polynomial order of the
%            approximation being performed. If the order is the same in all
%            dimensions, then a single scalar can be passed.
% numCentralKnots This optional parameter is used if one passes a higher
%            density of points in tau than is needed, because one wants to
%            perform least squares interpolation. This is the numDimsX1 or
%            1XnumDims list of the number of knots in the central region of
%            knots in each dimension. It varies from 0 to n-k. The higher
%            the number, the more knots are used total. A total of
%            numCentralKnots+2*k will be returned by this function, though,
%            as explained below, the first and last knots in each dimension
%            are not used.
%  useSparse The linear system being solved for the coefficients can be
%            quite large, but is not necessarily dense. If useSparse is
%            true, then the matrix for the system will be allocated as a
%            sparse matrix. This can slow down the algorithm, but it can
%            also keep Matlab from running out of memory on large systems.
%            The default if this parameter is omitted or an empty matrix is
%            passed is true.
%
%OUTPUTS: a A hypermatrix with numDim+1 indices containing the set of
%           coefficients for the b-splines covering all of the
%           interpolation intervals with the final index selecting the set
%           if numSets>1. If numSets=1, then the final index is unitary,
%           meaning that there are effectively only numDim dimensions. This
%           and the following outputs can be passed to
%           BSplineInterpValMultiDim to perform interpolation over the
%           region spanned by the tau values.
%         t The maxNumKnotsXnumDims set of knots for the interpolation
%           function. The first and last k(curDim)-1 knots in each
%           dimension are outside of the ends of the region with data or
%           mark the ends of the region with data. The number of rows is
%           the maximum number needed for all of the dimensions. tLength
%           says how many items in each column are actually used. The knots
%           are the same for each set.
%   tLength A numDimX1 vector where tLength(i) says the number of elements
%           in t(:,i). The lengths are the same for each set.
%
%The comments ot the function BSplinePolyFit describe how Equation 7 in
%Chapter IX is used in one dimension. The issue of multiple dimensions is
%addressed in Chapter 17 of [1] with a focus on solving the problem in
%2D. The basic equation to be solved in 2D is given by Equation 11 in
%Chapter 17 of [1] and it can be easily extended to an arbitrary number
%of dimensions. Chapter 17 gives a solution that is more efficient in 2D
%than solving Equation 11. However, we choose to solve the equation to make
%the problem valid for an arbitrary number of dimensions.
%
%The ability to do a least-squares fit involves making fewer inner knots
%than suggested in the text. Rather than making n-k knots, we make
%numCentralKnots. However, this invalidates the approximation in Equation
%10 of Chapter 13 of [1] for selecting knots. The correction applied is to
%make the average taken in the equation to involve more points.
%Additionally, the knots in Chapter 13 are univariate. We make them
%multivariate by taking a tensor product.
%
%EXAMPLE 1:
%Here, we fit a curve to 2D data on an irregular grid. We choose a
%sufficiently high order that it recreates the polynomial exactly and we
%verify this by looking at the maximum error on a fine grid and plotting
%the results.
% f=@(x,y)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y);
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
% 
% numPointsX=100;
% numPointsY=101;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% [X,Y]=ndgrid(pointsX,pointsY);
% x=[X(:).';Y(:).'];
% 
% zTrue=f(X,Y);
% 
% z=BSplineInterpValMultiDim(x,t,tLength,a,k);
% z=reshape(z,[numPointsX,numPointsY]);
% 
% figure(1)
% clf
% surface(X,Y,z,'EdgeColor','none')
% title('Interpolated Surface')
% 
% figure(2)
% clf
% surface(X,Y,zTrue,'EdgeColor','none')
% title('True Surface')
% max(abs(zTrue(:)-z(:)))
%One sees that the curves are the same and that the difference at any point
%is on the order of finite precision error. This is consistent with the
%fact that the order of the spline fit was chosen to be the same as the
%order of the true polynomial.
%
%EXAMPLE 2:
%Here, we fit a curve to 3D data on an irregular grid. We choose a
%sufficiently high order that it recreates the polynomial exactly and we
%verify this by looking at the maximum error on a fine grid.
% f=@(x,y,z)(x.^4-2*x.^2+x).*(y.^4-2*y.^2+y).*(z.^4-3*z.^2+z);
% numDims=3;
% numPointsX=10;
% numPointsY=11;
% numPointsZ=9;
% tauLengths=[numPointsX;numPointsY;numPointsZ];
% tau=NaN(max(tauLengths),numDims);
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
% 
% numPointsX=50;
% numPointsY=51;
% numPointsZ=49;
% pointsX=linspace(-1.5,1.5,numPointsX);
% pointsY=linspace(-1.5,1.5,numPointsY);
% pointsZ=linspace(-1.5,1.5,numPointsZ);
% [X,Y,Z]=ndgrid(pointsX,pointsY,pointsZ);
% x=[X(:).';Y(:).';Z(:).'];
% 
% wTrue=f(X,Y,Z);
% w=BSplineInterpValMultiDim(x,t,tLength,a,k);
% max(abs(wTrue(:)-w(:)))
%As was the case with Example 1, the maximum error is on the order of what
%one would expect due to finite precision limitations.
%
%EXAMPLE 3:
%In this case, a complex 2D function is fit and interpolated.
% numPoints=50;
% xMin=0;
% xMax=3;
% yMin=0;
% yMax=3;
% x=linspace(xMin,xMax,numPoints);
% y=linspace(yMin,yMax,numPoints);
% [X,Y]=ndgrid(x,y);
% f=@(x,y)((sin(5*y)+1j*cos(10*x)).*exp(1j*(x.^2-2*y.^2)));
% 
% fxy=f(X(:),Y(:));
% k=[5;5];
% tau=[x(:),y(:)];
% tauLengths=[numPoints,numPoints];
% [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,fxy,k);
% 
% %Interpolate on a finer grid.
% numPoints=100;
% x=linspace(xMin,xMax,numPoints);
% y=linspace(xMin,xMax,numPoints);
% [X,Y]=ndgrid(x,y);
% 
% fxy=f(X(:),Y(:));
% pts=[X(:).';Y(:).'];
% fxyInterp=reshape(BSplineInterpValMultiDim(pts,t,tLength,a,k),[numPoints,numPoints]);
% 
% fxy=reshape(fxy,[numPoints,numPoints]);
% 
% figure(1)
% clf
% hold on
% surface(X,Y,real(fxy),'EdgeColor','none')
% title('Real Function Values')
% colorbar()
% 
% figure(2)
% clf
% hold on
% surface(X,Y,imag(fxy),'EdgeColor','none')
% title('Imaginary Function Values')
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

numDims=size(tau,2);
n=tauLengths(:);
numSets=size(y,2);

if(nargin<6)
    useSparse=false;
end

%The total number of points.
NTotal=prod(n);
if(NTotal~=size(y,1))
    error('The dimensionalities of the inputs are not consistent.')
end

if(isscalar(k))
    k=repmat(k,[numDims,1]); 
end
k=k(:);

if(nargin<5||isempty(numCentralKnots))
    numCentralKnots=n-k;
elseif(any(numCentralKnots<0)||any(numCentralKnots(:)>n-k))
    error('The value of numCentralKnots is invalid.')
end
numCentralKnots=numCentralKnots(:);

deltaVal=n-k-numCentralKnots;

tLength=numCentralKnots+2*k;
maxTLength=max(tLength);

t=zeros(maxTLength,numDims);
for curDim=1:numDims
    taLengthCur=tauLengths(curDim);
    tauCur=tau(1:taLengthCur,curDim);
    kCur=k(curDim);
    deltaValCur=deltaVal(curDim);
    numCentKnotsCur=numCentralKnots(curDim);
    
    t(1:kCur,curDim)=tauCur(1);
    t((numCentKnotsCur+kCur+1):(numCentKnotsCur+2*kCur),curDim)=tauCur(n(curDim));

    %The ad-hoc point choice from 10 of Chapter 13, modified to take more
    %points if a least-squares solution is being performed.
    for i=1:numCentKnotsCur
        t(k+i,curDim)=sum(tauCur((i+1):(i+kCur-1+deltaValCur)))/(kCur-1+deltaValCur);
    end
end

%Given the locations of the knots, we build the equations. The a values are
%of the form a(d1,d2,..dnumDims). To access them, we stack them one after
%another in the same order as a(:) would in Matlab.
totalNumA=prod(numCentralKnots+k);

if(useSparse)
    %The amount of memory used will grow as B is filled in.
    B=sparse(NTotal,totalNumA);
else
    B=zeros(NTotal,totalNumA);
end

spanMin=ones(numDims,1);
tIdxVals=k;

%tauIdxVals holds the index of the value of tau used in each dimension
tauIdxVals=ones(numDims,1);
tauVal=zeros(numDims,1);
%Fill in tauVal with the current value.
for curDim=1:numDims
    tauVal(curDim)=tau(1,curDim);
end

%The number of B values produced per call to evalBSplinePolysMultiDim.
numBCurEls=prod(k);

%We have to go through all of the indices.
for curEq=1:NTotal
    BCur=evalBSplinePolysMultiDim(t,tLength,tauVal,k,tIdxVals);
    %We have to put the elements in BCur into the proper spots in B. These
    %spots start at spanMin in each dimension.
    for curEl=1:numBCurEls
        dimIdx=index2NDim(k,curEl);
        %The multidimensional index of the a value by which the current B
        %value is multiplied.
        aDimIdx=spanMin+dimIdx-1;
        aIdx=nDim2Index(n-deltaVal,aDimIdx);        
        
        B(curEq,aIdx)=BCur(curEl);
    end
    
    if(curEq==NTotal)
        break;
    end
    
    %Now, we increment tauIdxVals to move onto a new value of tau in each
    %dimension. the values are kept in the tauVal array. Each time a value
    %is incremented, we have to check whether the corresponding entry in
    %tIdxVals and spanMin should be incremented.
    curIdx=1;
    tauIdxVals(curIdx)=tauIdxVals(curIdx)+1;
    while(tauIdxVals(curIdx)>tauLengths(curIdx))
        %If it exceeds k, then go back to 1.
        tauIdxVals(curIdx)=1;
        tauVal(curIdx)=tau(1,curIdx);

        tIdxVals(curIdx)=k(curIdx);
        spanMin(curIdx)=1;
        
        %Increment tIdxVals and spanMin as necessary to reach the valid
        %region.
        while(tauVal(curIdx)>t(tIdxVals(curIdx)+1,curIdx))
            tIdxVals(curIdx)=tIdxVals(curIdx)+1;
            spanMin(curIdx)=spanMin(curIdx)+1;
        end
 
        curIdx=curIdx+1;
        tauIdxVals(curIdx)=tauIdxVals(curIdx)+1;
    end

    %Increment tIdxVals and spanMin as necessary to reach the valid
    %region.
    tauVal(curIdx)=tau(tauIdxVals(curIdx),curIdx);
    while(tauVal(curIdx)>t(tIdxVals(curIdx)+1,curIdx))
        tIdxVals(curIdx)=tIdxVals(curIdx)+1;
        spanMin(curIdx)=spanMin(curIdx)+1;
    end
end

a=B\y;
if(numDims>1)
    a=reshape(a,[(numCentralKnots+k)',numSets]);
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
