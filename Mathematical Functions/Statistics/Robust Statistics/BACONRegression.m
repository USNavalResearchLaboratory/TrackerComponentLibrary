function [betam,sigmam]=BACONRegression(x,y,m,hasIntercept,BACONParams)
%%BACONREGRESSION Perform multilinear regression with outliers by using the
%            blocked adaptive computationally effcient outlier nominators
%            (BACON) algorithm of [1]. This algorithm tries to detect and
%            eliminate outliers. The algorithm finds betam such that
%            y=betam'*x for an input x and y. If hasIntercept is true, then
%            the equation is y=betam'*[x;1] as the last element of betam
%            will be an additive constant.
%
%INPUTS: x A pXn set of n p-dimensional samples of the p-dimensional
%          independent parameter. It is required that n>p.
%        y A 1Xn or nX1 vector of the dependant parameter.
%        m The number of points to use in the ultimate set of points that
%          are used for regression. The default if omitted or an empty
%          matrix is passed is m=max(p+1,ceil(n*3/4));
% hasIntercept A boolean value indicating whether the fit has a nonzero y
%          intercept. The default if omitted or an empty matrix is passed
%          is true.
% BACONParameters The function meanCovBACON is used in the first step of
%          the regression algorithm. This structure can contain elements
%          'version', 'gammaVal', and 'maxIter' to override the defaults of
%          the meanCovBACON function.
%
%OUTPUTS: betam The regression coefficient as described above.
%        sigma2 The residual mean.  This is SSE/(n-p) (where p is actually
%               one larger than the number of rows of x if
%               hasIntercept=true) and where SSE is the vector of ordinary
%               residuals e'*e with e=y(:)-x'*betam. x is augmented by a
%               row of ones if a nonzero y-intercept is desired.
%
%This function implements Algorithms 4 and 5 of [1].
%
%EXAMPLE:
%Here, we show that given noise-free samples of a line polluted 10% by
%samples from another line, one is able to find the dominant line exactly.
% numTrue=100;
% betaTrue=[14;-8];
% yIntTrue=42;
% xTrue=rand(2,numTrue);
% yTrue=sum(bsxfun(@times,betaTrue,xTrue),1)+yIntTrue;
% numBad=10;
% betaBad=[-80;32];
% yIntBad=pi;
% xBad=rand(2,numBad);
% yBad=sum(bsxfun(@times,betaBad,xBad),1)+yIntBad;
% 
% x=[xTrue,xBad];
% y=[yTrue,yBad];
% numTotal=size(x,2);
% 
% m=fix(0.8*numTotal);
% [betam,sigmam]=BACONRegression(x,y,m,true)
%One will find that betam=[14;-8] within expected finite precision
%limitations. Thus, the line is fit despite the outliers.
%
%REFERENCES:
%[1] N. Billor, A. S. Hadi, and P. F. Velleman, "BACON: Blocked adaptive
%    computationally efficient outlier nominators," Computational
%    Statistics and Data Analysis, vol. 34, no. 3, pp. 279-298, 28 Sep.
%    2000.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

p=size(x,1);
n=size(x,2);

maxIter=[];
version=[];
gammaVal=[];
if(nargin>4&&~isempty(BACONParams))
    if(isfield(BACONParams,'maxIter'))
        maxIter=BACONParams.maxIter;
    end
    
    if(isfield(BACONParams,'version'))
        version=BACONParams.version;
    end
    
    if(isfield(BACONParams,'gammaVal'))
        gammaVal=BACONParams.gammaVal;
    end
end

if(nargin<4||isempty(hasIntercept))
    hasIntercept=true; 
end

if(nargin<3||isempty(m))
    m=max(p+1,ceil(n*3/4));
end

y=y(:);

if(rank(x)<p)
    error('x must be rank p.');
end

%Allocate the maximum space to hold the subset matrix and the indices.
Xm=zeros(p,n);
idxVals=zeros(1,n);

%Algorithm 4, Step 0
[~,~,d,idx]=meanCovBACON(x,m,version,gammaVal,maxIter);

subsetSize=length(idx);
idxVals(1:subsetSize)=idx;

Xm(1:p,1:subsetSize)=x(:,idxVals(1:subsetSize));
if(rank(Xm(1:p,1:subsetSize))<p&&subsetSize<n)
    [Xm,subsetSize,idxVals]=addObs2MakeFullRank(x,d,Xm,subsetSize,idxVals);
end

if(hasIntercept)
    Xm(p+1,:)=1;%Add another row.
    x(p+1,:)=1;
    p=p+1;
end

[tAbs,betam,sigmam]=regressionValsAndCost(x,y,Xm(:,1:subsetSize),idxVals(1:subsetSize));

subsetSize=p+1;
[~,idx]=mink(tAbs,subsetSize);
idxVals(1:subsetSize)=idx;
Xm(:,1:subsetSize)=x(:,idxVals(1:subsetSize));
while(subsetSize<m)
    %Algorithm 4 Step 1
    if(rank(Xm(:,1:subsetSize))<p&&subsetSize<n)
        [Xm,subsetSizeNew,idxVals]=addObs2MakeFullRank(x,tAbs,Xm,subsetSize,idxVals);
        m=min(n,m+(subsetSizeNew-subsetSize));
        subsetSize=subsetSizeNew;
    end

    [tAbs,betam,sigmam]=regressionValsAndCost(x,y,Xm(:,1:subsetSize),idxVals(1:subsetSize));

    subsetSize=min(n,subsetSize+1);
    [~,idx]=mink(tAbs,subsetSize);
    idxVals(1:subsetSize)=idx;
    Xm(:,1:subsetSize)=x(:,idxVals(1:subsetSize));
end
end

function [Xm,subsetSize,idxVals]=addObs2MakeFullRank(x,d,Xm,subsetSize,idxVals)
%%ADDOBS2MAKEFULLRANK If a given subset of observations is not full rank,
%           this function keeps adding observations until it is full rank.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

p=size(x,1);
n=size(x,2);

sel=ones(1,n);
sel(idxVals(1:subsetSize))=false;
idxOmit=find(sel);

dOmit=d(idxOmit);
[~,idx]=sort(dOmit,'ascend');
idxOmit=idxOmit(idx);

cur2Add=1;
while(subsetSize<n)
    subsetSize=subsetSize+1;
    idxVals(subsetSize)=idxOmit(cur2Add);
    Xm(:,subsetSize)=x(idxOmit(cur2Add));

    %If it has become observable.
    if(~rank(Xm(:,1:subsetSize))<p)
        break;
    end

    cur2Add=cur2Add+1;
end
end

function [tAbs,betam,sigmam]=regressionValsAndCost(x,y,Xm,idxVals)
%%REGRESSIONVALSANDCOST This function is used to essentially fit the line
%           and evaluate Equations 5 and 6 in algorithms 4 and 5  of [1].
%
%REFERENCES:
%[1] N. Billor, A. S. Hadi, and P. F. Velleman, "BACON: Blocked adaptive
%    computationally efficient outlier nominators," Computational
%    Statistics and Data Analysis, vol. 34, no. 3, pp. 279-298, 28 Sep.
%    2000.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=size(x,2);

[betam,sigmam2]=multilinRegress(Xm,y(idxVals),[],false);
sigmam=sqrt(sigmam2);

selSet=false(n,1);
selSet(idxVals)=true;

RxInv=inv(Xm*Xm');

tAbs=zeros(n,1);
for i=1:n
    %Equation 5
    if(selSet(i))
        tAbs(i)=abs(y(i)-x(:,i)'*betam)/(sigmam*sqrt(abs(1-x(:,i)'*RxInv*x(:,i))));
    else
        tAbs(i)=abs(y(i)-x(:,i)'*betam)/(sigmam*sqrt(1+x(:,i)'*RxInv*x(:,i)));
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
