function [CMerged,omega,xMerged]=covarianceIntersect(C1,C2,optCrit,x1,x2)
%%COVARIANCEINTERSECT Perform covariance intersection. This is a method of
%                     fusing the first two moments of estimates when the
%                     correlation between the estimates is unknown. Given
%                     the covariance matrices of the estimates, this
%                     function returns a covariance matrix and a scaling
%                     factor, which can be used to fuse estimates.
%                     Alternatively, if the estimates themselves are given,
%                     this function can also return the merged estimate.
%                     Compare this function to ellipsoidIntersect.
%
%INPUTS: C1 The nXn positive definite covariance matrix of the first
%           estimate.
%        C2 The nXn positive definite covariance matrix of the second
%           estimate.
%   optCrit An optional  parameter specifying the optimization criterion 
%           for the fusion. Possible values are
%           'det' (The default if omitted or an empty matrix is passed) The
%                 fusion is performed to optimize the determinant of the
%                 matrix on the output.
%           'trace' The fusion is performed to minimize the trace of the
%                 output matrix.
%    x1,x2 The optional nX1 estimate vectors to be merged. These are only
%          used if xMerged is requested on the output.
%
%OUTPUTS: CMerged The nXn fused covariance matrix.
%           omega The weighting that played a role in the fusion of the
%                 covariance matrix/ state estimates. Specifically,
%                 CMerged=inv(omega*inv(C1)+(1-omega)*inv(C2))
%                 and the merged estimates are
%                 xMerged=CMerged*(omega*inv(C1)*x1+(1-omega)*inv(C2)*x2)
%         xMerged The merged estimate. This requires x1,x2 to be given on
%                 the input. The covariance matrix associated with the
%                 merged estimate is CMerged.
%
%The first (non-dissertation) publication of covariance intersection for
%fusing measurements having unknown correlations is [2]. Various optimality
%criteria are derived in [3]. However, the most practical approach is that
%of [1], which clearly expresses the solution in terms of the solution of
%polynomials, with some solutions given in closed form. Thus, this function
%uses the closed form solutions, when available, and the general polynomial
%solutions otherwise.
%
%Note that covariance intersection is overly convervative in its covariance
%estimates as mentioned in Chapter 9.3.7 of [5].
%
%EXAMPLE:
%Here we use the numerical values form the example for ellipsoidal
%intersection in [1]. The plots in the paper appear to be incorrect.
% xi=[1;-2];
% Pi=[3,0;0,0.4];
% xj=[-2;-1];
% Pj=[2,-0.8;-0.8,1];
% 
% [P,~,x]=covarianceIntersect(Pi,Pj,'det',xi,xj);
% figure()
% clf
% hold on
% drawEllipse(xi,inv(Pi),[],'--r')
% drawEllipse(xj,inv(Pj),[],'--g')
% drawEllipse(x,inv(P),[],'-b')
%
%REFERENCES:
%[1] M. Reinhardt, B. Noack, and U. D. Hanebeck, "Closed-form optimization
%    of covariance intersection for low-dimensional matrices," in 
%    Proceedings of the 15th International Conference on Information
%    Fusion, Singapore, 9-12 Jun. 2012, pp. 1891-1896.
%[2] S. J. Julier and J. K. Uhlmann, "A non-divergent estimation algorithm
%    in the presence of unknown correlations," in Proceedings of the
%    American Control Conference, vol. 4, Albuquerque, NM, 4-6 Jun. 1997,
%    pp. 2369-2373.
%[3] L. Chen, P. O. Arambel, and R. K. Mehra, "Estimation under unknown
%    correlation: Covariance intersection revisited," IEEE Transactions on
%    Automatic Control, vol. 47, no. 11, pp. 1879-1882, Nov. 2002.
%[4] J. Sijs, M. Lazar, and P. Bosch, "State fusion with unknown
%    correlation: Ellipsoidal intersection," in Proceedings of the 2010
%    American Control Conference, Baltimore, MD, 30 Jun. - 2 Jul. 2010.
%[5] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(optCrit))
    optCrit='det';
end

n=length(C1);

%The transformation matrix in Equation 2 of [1], but here we use an
%SVD-based algorithm instead of an eigenvalue-based algorithm, because it
%is numerically more stable.
A=twoMatDiag(C1,C2);

%The positive scalars d.
d=diag(A*C2*A');
dBar=1./d;
dTilde=dBar./(1-dBar);

invC1=inv(C1);
invC2=inv(C2);

%Determinant minimization
switch(optCrit)
    case 'det'
        if(n==1)
            if(C1<C2)
                omega=1;
            else
                omega=0;
            end
        elseif(n==2)
            %The solution is Equation 11 in [1].
            omega=-(1/2)*(dTilde(1)+dTilde(2));
        elseif(n==3)
            sqrtTerm=sqrt(sum(dTilde.^2)-dTilde(1)*dTilde(2)-dTilde(1)*dTilde(3)-dTilde(2)*dTilde(3));
            %Equation 12 in [1].
            omega(1,1)=-(1/3)*(sum(dTilde)+sqrtTerm);
            omega(2,1)=-(1/3)*(sum(dTilde)-sqrtTerm);
        elseif(n==4)
            b=(3/4)*sum(dTilde);

            c=0;
            for i=1:n
               for j=(i+1):n
                   c=c+dTilde(i)*dTilde(j);
               end
            end
            c=(1/2)*c;

            d=0;
            for i=1:n
                for j=(i+1):n
                    for k=(j+1):n
                        d=d+dTilde(i)*dTilde(j)*dTilde(k);
                    end
                end
            end
            d=(1/4)*d;

            Q=sqrt((2*b^3-9*b*c+27*d)^2-4*(b^2-3*c)^3);
            C=((1/2)*(Q+2*b^3-9*b*c+27*d))^(1/3);

            %Equation 20
            omega(1,1)=-(b/3)-C/3-(b^2-3*c)/(3*C);
            omega(2,1)=-(b/3)+C*(1+1j*sqrt(3))/6+(1-1j*sqrt(3))*(b^2-3*c)/(6*C);
            omega(3,1)=-(b/3)+C*(1-1j*sqrt(3))/6+(1+1j*sqrt(3))*(b^2-3*c)/(6*C);
        else%For n>4
            %The pi values are given in the unnumbered equation after
            %Equation 10 in [1].
            piVals=zeros(1,n-1);
            for curPi=1:(n-1)
                piIdx=1:curPi;
                while(1)
                    %Add the current  product
                    piVals(curPi)=piVals(curPi)+prod(dTilde(piIdx));
                    
                    %Increment the value in the innermost sum.
                    curLevel=curPi;
                    allDone=false;
                    while(piIdx(curLevel)+1>(n-curPi+curLevel))
                        curLevel=curLevel-1;
                        if(curLevel<1)
                            allDone=true;
                            break;
                        end
                    end
                    if(allDone)
                        break;
                    end
                    %Increment the index
                    piIdx(curLevel)=piIdx(curLevel)+1;
                    
                    %Fill in the minimum values in further nested sums.
                    for k=(curLevel+1):curPi
                        piIdx(k)=piIdx(k-1)+1;
                    end
                end
            end
            
            %Given the pi values, build the polynomial in omega.
            poly2Solve=(n:-1:1).*[1,piVals];
            
            %Omega hypotheses are the roots of an order n-1 polynomial in
            %Equation 10 of [1].
            omega=roots(poly2Solve);
        end
        
        if(n>1)
        %Given a set of omega values, we have to choose the best solution.
        %We also have to include omega=0 and omega=1 as possibilities as
        %per Algorithm 1 in [1]. Some of the solutions might be slightly
        %imaginary due to finite precision errors. We will just discard the
        %imaginary parts. This does not affect the choice of the optimal
        %solution.
            omega=[real(omega);0;1];
            numHyp=length(omega);

            minCost=Inf;
            minOmega=[];
            for curHyp=1:numHyp
                if(omega(curHyp)>=0&&omega(curHyp)<=1)
                    %The determinant of a matrix inverse is the inverse of
                    %the determinant.
                    cost=1/det(omega(curHyp)*invC1+(1-omega(curHyp))*invC2);

                    if(cost<minCost)
                       minCost=cost;
                       minOmega=omega(curHyp);
                    end
                end
            end

            omega=minOmega;
        end
    case 'trace'
        if(n==1)
            %The solution is the same as for the determinant minimization
            %criterion.
            if(C1<C2)
                omega=1;
            else
                omega=0;
            end
        elseif(n==2)
            invA=inv(A);
            a=diag(invA'*invA);
            
            denom=a(1)*(1+dTilde(1))+a(2)*(1+dTilde(2));
            p=(a(1)*dTilde(2)*(1+dTilde(1))+a(2)*dTilde(1)*(1+dTilde(2)))/denom;
            q=(a(1)*dTilde(2)^2*(1+dTilde(1))+a(2)*dTilde(1)^2*(1+dTilde(2)))/denom;
            
            %Equation 19
            omega(1,1)=-p+sqrt(p^2-q);
            omega(2,1)=-p-sqrt(p^2-q);
        else%for n>=3
            %Use convolutions to build up the polynomial in Equation 17 in
            %[1].
            invA=inv(A);
            a=diag(invA'*invA);
            
            %The polynomial is degree 2*(n-1).
            numPolyEls=2*(n-1)+1;
            poly2Solve=zeros(numPolyEls,1);
            for i=1:n
                coeff=a(i)*(1+dTilde(i));
                
                jPoly=zeros(numPolyEls,1);
                jPoly(end)=1;
                for j=1:n
                    if(j==i)
                        continue;
                    end
                    
                    jPoly=conv(jPoly,[1;2*dTilde(j);dTilde(j)^2]);
                    %Get rid of zero-padding at the beginning.
                    jPoly=jPoly((end-numPolyEls+1):end);
                end
                poly2Solve=poly2Solve+coeff*jPoly;
            end
            
            %omega hypotheses are the roots of the equation.
            omega=roots(poly2Solve);
        end
        
        %Given a set of omega values, we have to choose the best solution.
        %We also have to include omega=0 and omega=1 as possibilities as
        %per Algorithm 1 in [1]. Some of the solutions might be slightly
        %imaginary due to finite precision errors. We will just discard the
        %imaginary parts. This does not affect the choice of the optimal
        %solution.
        if(n>1)
            omega=[real(omega);0;1];
            numHyp=length(omega);

            minCost=Inf;
            minOmega=[];
            for curHyp=1:numHyp
                if(omega(curHyp)>=0&&omega(curHyp)<=1)
                    cost=trace(inv(omega(curHyp)*invC1+(1-omega(curHyp))*invC2));

                    if(cost<minCost)
                       minCost=cost;
                       minOmega=omega(curHyp);
                    end
                end
            end

            omega=minOmega;
        end
    otherwise
        error('Invalid optimality criterion chosen.')
end

CMergedInv=omega*invC1+(1-omega)*invC2;
CMerged=inv(CMergedInv);

if(nargout>2)
    xMerged=CMergedInv\(omega*invC1*x1+(1-omega)*invC2*x2);
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
