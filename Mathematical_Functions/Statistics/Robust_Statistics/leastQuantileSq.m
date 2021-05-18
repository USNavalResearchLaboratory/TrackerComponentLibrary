function [thetaCMin,b,r2qMin]=leastQuantileSq(x,y,q,hasIntercept)
%%LEASTQUANTILESQ Compute the least quantile of squares estimate (a robust
%     type of multilinear regression). The least median of squares
%     estimate is a special case with q=ceil(n/2)+ceil((xDim+1)/2). Given a
%     set of n xDimX1 vectors x and n scalar values of y, this function
%     finds the xDimX1 vector theta such that
%     thetac=minimize_theta ( median_i (y_i-x_i'*theta+b )^2 )
%     whereby for median the qth order statistic is used. 
%
%INPUTS: x The xDimXn set of n values of the independent parameter. Note
%          that (p+1)>=n.
%        y The 1Xn or nX1 set of values of the (scalar) dependent
%          parameter.
%        q The (integer index) value of the sample to use for the median.
%          If omitted or an empty matrix is passed, then
%          q=ceil(n/2)+ceil((xDim+1)/2) is used, which corresponds to the
%          least median of squares estimate.
% hasIntercept If true, then the estimation problem is assumed to have a
%           nonzero y intercept, which should be estimated. Otherwise, b is
%           returned as an empty matrix. The default if omitted or an empty
%           matrix is passed is false.
%
%OUTPUTS: thetaCMin The xDimX1 slope vector.
%                 b If hasIntercept=true, this is the y intercept of the
%                   estimate. otherwise, this is an empty matrix.
%            r2qMin The median value of the squared residuals. In other
%                   words, the value of the cost function given thetaCMin
%                   and b.
%
%The least median of squares regression is introduced in [1], though no
%specific implementation details are given. The implementation used here,
%is based on what is in [2]. The algorithm assumes that every set of xDim+1
%columns of x is full row rank (rank xDim). If not, then it is not
%guaranteed that the solution found is globaly optimal. The computational
%complexity of the algorithm scales as binomial(n,xDim+1).
%
%EXAMPLE:
%This example uses the salinity (as x) and H20 Flow (as y) from the clouse
%seeding data that is given in [3]. We plot the lines.
% x=[7.6;7.7;4.3;5.9;5.0;6.5;8.3;8.2;13.2;12.6;10.4;10.8;13.1;12.3;10.4;10.5;7.7;9.5;12.0;12.6;13.6;14.1;13.5;11.5;12.0;13.0;14.1;15.1]';
% y=[23.005;23.873;26.417;24.868;29.895;24.200;23.215;21.862;22.274;23.830;25.144;22.430;21.785;22.380;23.927;33.443;24.859;22.686;21.789;22.041;21.033;21.005;25.865;26.290;22.932;21.313;20.769;21.393];
% 
% hasIntercept=true;
% [thetaCMin,b]=leastQuantileSq(x,y,[],hasIntercept);
% coeffs=multilinRegress([x;y']);
% 
% figure(1)
% clf
% hold on
% scatter(x,y,'filled')
% numPoints=100;
% xTest=linspace(4,16,numPoints);
% yLME=xTest'*thetaCMin+b;
% yLS=xTest'*coeffs(1:(end-1))+coeffs(end);
% plot(xTest,yLME,'-b','linewidth',2);
% plot(xTest,yLS,'-k','linewidth',2);
% h1=xlabel('Salinity');
% h2=ylabel('H20 Flow');
% legend('Data Point','Least Median of Squares','Least Squares','location','NorthEast')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] P. J. Rousseeuw, "Least median of squares regression," Journal of the
%    American Statistical Association, vol. 79, no. 388, pp. 871-880, Dec.
%    1984.
%[2] A. J. Stromberg, "Computing the exact least median of squares estimate
%    and stability diagnostics in multiple linear regression," SIAM Journal
%    on Scientific Computing, vol. 14, no. 6, pp. 1289-1299, Nov. 1993.
%[3] D. Ruppert and R. J. Carroll, "Trimmed least squares estimation in the
%    linear model," Journal of the American Statistical Association, vol.
%    75, no. 372, pp. 828-838, Dec. 1980.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

p=size(x,1);%Dimensionality of vectors.
n=size(x,2);%Number of samples.

if(n<p+1)
    error('It is required that n>=p+1');
end

if(nargin<4||isempty(hasIntercept))
    hasIntercept=false;
end

if(hasIntercept)
   p=p+1;
   x=[x;ones(1,n)];
end

%The default q is for the least median of squares estimate.
if(nargin<3||isempty(q))
    q=ceil(n/2)+ceil((p+1)/2);
end
numQ=length(q);
y=y(:);

r2qMin=Inf(numQ,1);
thetaCMin=NaN(p,numQ);

%Go through all combinations of p+1 points.
curCombo=0:p;
while(~isempty(curCombo))
    comboIdx=curCombo+1;
    
    XCur=x(:,comboIdx);
    yCur=y(comboIdx);
    
    %The method from Theorem 2:
    M=pinv(XCur');
    thetaLS=M*yCur;
    
    %Get the LS residuals
    r=yCur-XCur'*thetaLS;

    %Adjust the estimate.
    epsilon=sum(r.^2)/sum(abs(r));
    s=sign(r);
    thetaC=M*(yCur-epsilon*s);
    
    %Get the residuals from the new estimate.
    r=y-x'*thetaC;
    r2=sort(r.^2,'ascend');
    
    for curQ=1:numQ
        if(r2(q(curQ))<r2qMin(curQ))
            r2qMin(curQ)=r2(q);
            thetaCMin(:,curQ)=thetaC;
        end
    end

    curCombo=getNextCombo(curCombo,n);
end

if(hasIntercept)
    b=thetaCMin(end,:);
    thetaCMin=thetaCMin(1:(p-1),:);
else
    b=[];
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
