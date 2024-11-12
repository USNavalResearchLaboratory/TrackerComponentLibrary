function [Q,b,c,xMode]=secondOrderPolyFitExtremum(yVals,xVals,zeroOffDiags,transformX)
%%SECONDORDERPOLYFITEXTREMUM Given a function taking a vector input and
%   outputting a scalar value, this function fits a multivariate quadratic
%   equation of the form y=x'*Q*x+b'*x+c. The function outputs the expremum
%   of this fitted function, if desired. This can be useful for taking
%   points around a sampled peak and interpolating a more accurate peak.
%
%INPUTS: yVals A length numMeas vector of the real function values
%              evaluated at the points in xVals.
%        xVals A numDimXnumMeas matrix of the real vector points at which
%              the values in yVals were obtained.
% zeroOffDiags If this is true, then it is assumed that the off-diagonal
%              terms of Q are zero. The default if omitted or an empty
%              matrix is passed is false.
%   transformX If this is true, then the element of xVals are transformed
%              prior to fitting so that they are centered about the origin
%              and each dimension has the same maximum magnitude. After
%              fitting, the results are transformed back to the original
%              problem. This can be useful when performing least squares
%              fitting when the scale in the dimensions of x differ
%              notably. The default if omitted or an empty matrix is passed
%              is false.
%
%OUTPUTS: Q, b, c The numDimXnumDim, numDimX1 and 1X1 fitted parameters in
%                 the equation y=x'*Q*x+b'*x+c.
%           xMode The extremum of the fitted function.
%
%Equation pair of yVals(i) and xVals(:,i) defines a linear equation in
%terms of the unknown parameters of the quadratic equation. This function
%just solves that equation, using pinv to deal with being given more points
%than necessary.
%
%For there to be a unique solution, if transformX=false,
%3*binomial(2+numDim,numDim-1)/numDim measurments are needed. If
%transformX=true, then 2*numDim+1 measurements is are needed.
%
%EXAMPLE 1:
%In this example, we take an actual quadratic equation in 2D, generate
%points, and demonstrate that within reasonable finite precision limits, one
%gets back the original Q, b, and c values. Also, the function is plotted
%and it can be seen that the minimum coincides with the actual function
%minimum.
% Q=[4,-1;
%    -1,2];
% b=[4;9];
% c=12;
% %The minimum number of points needed is 6. We are taking 9.
% xVals=[0,1,0,-1, 0,1,-1, 1,-1;
%        0,0,1, 0,-1,1, 1,-1,-1];
% numVals=size(xVals,2);
% yVals=zeros(numVals,1);
% for k=1:numVals
%     yVals(k)=xVals(:,k)'*Q*xVals(:,k)+b'*xVals(:,k)+c;
% end
% [Q1,b1,c1,xMode1]=secondOrderPolyFitExtremum(yVals,xVals);
% v1=[Q(:);b(:);c];
% v2=[Q1(:);b1(:);c1];
% RelErr=max(abs((v2-v1)./v1))
% 
% numPts1=100;
% numPts2=101;
% xLin1=linspace(xMode1(1)-20,xMode1(1)+20,numPts1);
% xLin2=linspace(xMode1(2)-20,xMode1(2)+20,numPts2);
% [X1,X2]=meshgrid(xLin1,xLin2);
% x=[X1(:).';X2(:).'];
% numPts=numPts1*numPts2;
% Z=zeros(numPts2,numPts1);
% for curPt=1:numPts
%     xCur=x(:,curPt);
%     Z(curPt)=xCur'*Q*xCur+b'*xCur+c;
% end
% figure(1)
% clf
% hold on
% surface(X1,X2,Z,'EdgeColor','none')
% scatter(xMode1(1),xMode1(2),200,'.r')
% colormap(jet(512))
%
%EXAMPLE 2:
%This is the same as the first example without the plot, except the
%diagonals of the Q matrix are 0 and we are telling the function to assume
%that they are 0. In this instance, only 5 points are needed not 6, so only
%5 points are provided. The reconstructed parameters are within what would
%expect for finite precision limitations.
% Q=[4,0;
%    0,2];
% b=[4;9];
% c=12;
% %Taking the minimum number of points when assuming that the cross terms are
% %0, which is 5.
% xVals=[0,1,0,-1, 0
%        0,0,1, 0,-1];
% numVals=size(xVals,2);
% yVals=zeros(numVals,1);
% for k=1:numVals
%     yVals(k)=xVals(:,k)'*Q*xVals(:,k)+b'*xVals(:,k)+c;
% end
% zeroOffDiags=true;
% [Q1,b1,c1]=secondOrderPolyFitExtremum(yVals,xVals,zeroOffDiags);
% v1=[Q(:);b(:);c];
% v2=[Q1(:);b1(:);c1];
% RelErr=max(abs((v2-v1)./v1))
%
%April 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(transformX))
    transformX=false;
end

if(nargin<3||isempty(zeroOffDiags))
    zeroOffDiags=false;
end

numDim=size(xVals,1);
numMeas=length(yVals);

if(transformX)
    xBar=mean(xVals,2);
    xVals=bsxfun(@minus,xVals,xBar);
    xScal=max(abs(xVals),[],2);
    xScalInv=1./xScal;
    %This keeps the maximum magnitude value in each dimension=+/-1.
    xVals=bsxfun(@times,xVals,xScalInv);
end

if(zeroOffDiags)
    numCoeff=2*numDim+1;

    X=zeros(numMeas,numCoeff);
    for curMeas=1:numMeas
        x=xVals(:,curMeas);
        %The coefficients for the diagonal elements of Q.
        for i1=1:numDim
            X(curMeas,i1)=x(i1)^2;
        end
        %The coefficients for b.
        X(curMeas,(numDim+1):(2*numDim))=x;
        %The coefficient for c.
        X(curMeas,numCoeff)=1;
    end

    coeffs=pinv(X)*yVals(:);
    Q=diag(coeffs(1:numDim));
    b=coeffs((numDim+1):(2*numDim));
    c=coeffs(numCoeff);
else
    numCoeff=3*binomial(2+numDim,numDim-1)/numDim;
    
    X=zeros(numMeas,numCoeff);
    for curMeas=1:numMeas
        %First, the coefficients for the Q matrix.
        x=xVals(:,curMeas);
        idx=1;
        for i1=1:numDim
            %The Q(i,i) coefficient.
            X(curMeas,idx)=x(i1)^2;
            idx=idx+1;
            for i2=(i1+1):numDim
                %The Q(i,j) coefficient
                X(curMeas,idx)=2*x(i1)*x(i2);
    
                idx=idx+1;
            end
        end
        %Next, the coefficients for the b vector.
        for i1=1:numDim
            X(curMeas,idx)=x(i1);
            idx=idx+1;
        end
    
        %Finally, the c constant coefficient.
        X(curMeas,idx)=1;
    end
    
    coeffs=pinv(X)*yVals(:);
    
    %Break into the Q, b, and c values.
    Q=zeros(numDim,numDim);
    idx=1;
    for i1=1:numDim
        %The Q(i,i) coefficient.
        Q(i1,i1)=coeffs(idx);
        idx=idx+1;
        for i2=(i1+1):numDim
            %The Q(i,j) coefficient
            Q(i1,i2)=coeffs(idx);
            Q(i2,i1)=Q(i1,i2);
            idx=idx+1;
        end
    end
    
    b=coeffs(idx:(numCoeff-1));
    c=coeffs(numCoeff);
end

if(nargout>3)
    xMode=-(1/2)*(Q\b);
end

if(transformX)
    %Un-transform everything.
    if(nargout>3)
        xMode=xScal.*xMode+xBar;
    end

    Q=bsxfun(@times,bsxfun(@times,Q,xScalInv),xScalInv.');
    b=b.*xScalInv;
    c=c-b'*xBar+xBar'*Q*xBar;
    b=b-2*Q*xBar;
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
