function [theta,b,QMin]=LTSRegression(x,y,h)
%%LTSREGRESSION Perform least trimmed squares (LTS) regression using the
%               fast, approximate algorithm of [1] that is identified as
%               best for <600 data points. This means find theta and b that
%               minimizes 
%               sum_{i=1}^h r_i^2
%               where r_i is the ith largest value of
%               abs(y(i)-theta*x(:,i)-b)
%               This type of regression is robust to outliers.
%
%INPUTS: x,y These are the real samples of data. x is an xDimXN set of N
%            vectors and y is a set of corresponding scalar values. The
%            regression searches for theta and b to fit a line such that
%            y=theta'*x(:,i)+b.
%          h The trim value of the estimator. This must be between
%            fix((N+xDim+2)/2) and n. The default if omitted or an empty
%            matrix is passed is the value that allows for the highest
%            breakdown value, which is fix((N+xDim+2)/2).
%
%OUTPUTS: theta An xDimX1 vector. The multivariate slope of the line.
%             b The scalar y intercept of the data.
%          QMin The residual of the data.
%
%The algorithm described in [1] without the y intercept adjustment is used.
%
%EXAMPLE:
%This is similar to an example in [1]. There is a distribution that one
%might expect to be well approximated with a line corrupted by a minority
%of samples of a different distribution at the side. The plot is zoomed in
%so not all samples are shown. The standard linear regression result is
%shown for comparison.
% x=100*randn(1,800/2);
% y=x+1+randn(1,800/2);
% 
% xy=GaussianD.rand(200/2,[50;0],25*eye(2,2));
% x=[x,xy(1,:)];
% y=[y,xy(2,:)];
% 
% figure(1)
% clf
% hold on
% scatter(x,y,'b','filled')
% 
% coeffs=multilinRegress([x;y]);
% numPoints=100;
% xEst=linspace(-200,200,numPoints);
% yEst=xEst*coeffs(1)+coeffs(2);
% plot(xEst,yEst,'--k','linewidth',2)
% [theta,b]=LTSRegression(x,y);
% yEst=theta*xEst+b;
% plot(xEst,yEst,'-r','linewidth',2)
% axis([-30-20,70+20,-30-20,50+20])%Zoom in.
% legend('Samples','Least Squares Fit','LTS Fit','location','NorthWest')
% h1=xlabel('x');
% h2=ylabel('y');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] P. J. Rousseeuw and V. Katrien, "Computing LTS regression for large
%    data sets," Data Mining and Knowledge Discovery, vol. 12, no. 1, pp.
%    26-45, Jan. 2006.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(x,2);

x=[x;ones(1,n)];

y=y(:);

p=size(x,1);%p>=2

if(nargin<3||isempty(h))
    h=fix((n+p+1)/2);
end

if(h<fix((n+p+1)/2)||h>n)
   error('Invalid value of h provided.') 
end

totalStarts=binomial(n,p);
Q3Heap=BinaryHeap(10,false);
for k=1:min(totalStarts,500)
    if(totalStarts<500)
        %Do all possible initializations.
        idx=unrankCombination(k-1,n,p)+1;
    else        
        %Get a random subset of points in x.
        idx=randCombination(n,p)+1;
    end

    xCur=x(:,idx);
    yCur=y(idx);
    %If what we get is not full rank, then we just use a
    %pseudoinverse. This skips issues of determining which and how
    %many additional columns to add.
    if(rank(xCur)<p)
        theta0=pinv(xCur)'*yCur;
    else
        theta0=xCur'\yCur;
    end

    for i=1:3
        r0=y'-sum(bsxfun(@times,x,theta0),1);
        [~,idx]=sort(abs(r0),'ascend');
        idx=idx(1:h);

        xCur=x(:,idx);
        yCur=y(idx);
        theta0=pinv(xCur)'*yCur;
    end
    
    r2=(y'-sum(bsxfun(@times,x,theta0),1)).^2;
    Q3=sum(mink(r2,h));

    if(Q3Heap.heapSize()<10)
        Q3Heap.insert(Q3,idx);
    else
        topEl=Q3Heap.getTop();
        if(topEl.key>Q3)
            Q3Heap.deleteTop();
            Q3Heap.insert(Q3,idx);
        end
    end
end

%For the 10 Q3 values in the heap, iterate the C-step in Section 2
%until convergence.
theta0List=zeros(p,10);
QList=zeros(10,1);
for k=1:10
    topEl=Q3Heap.getTop();
    Q3Heap.deleteTop();

    idx=topEl.value;
    xCur=x(:,idx);
    yCur=y(idx);
    theta0=pinv(xCur)'*yCur;
    idxOld=zeros(h,1);

    %It was noted that convergence should be for fewer than 10
    %iterations, so we make the maximum 10.
    for curIter=1:10       
        if(all(idx(:)==idxOld(:)))
            break;
        end
        r0=y'-sum(bsxfun(@times,x,theta0),1);
        [~,idx]=sort(abs(r0),'ascend');
        idx=idx(1:h);

        xCur=x(:,idx);
        yCur=y(idx);
        theta0=pinv(xCur)'*yCur;
    end

    r2=(y'-sum(bsxfun(@times,x,theta0),1)).^2;
    Q=sum(mink(r2,h));
    
    QList(k)=Q;
    theta0List(:,k)=theta0;
end

[QMin,idx]=min(QList);
theta0Min=theta0List(:,idx);

theta=theta0Min(1:(p-1));
b=theta0Min(p);
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
