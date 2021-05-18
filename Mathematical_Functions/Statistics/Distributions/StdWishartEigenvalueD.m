classdef StdWishartEigenvalueD
%%STDWISHARTEIGENVALUED Functions to handle the distribution of the
%           eigenvalues of a real Wishart-distributed matrix with a scale
%           matrix consisting of the identity matrix multiplied by lambda
%           and a number of degrees of freedom of nu. This type of
%           distribution arises in the determination of confidence
%           intervals involved with extended target tracking in [1], and
%           [2].
%Implemented methods are: PDF, maxMinEigRegion2D, CDF2DMax
%
%REFERENCES:
%[1] M. Feldmann, D. Fränken, and W. Koch, "Tracking of extended objects
%    and group targets using random matrices," IEEE Transactions on Signal
%    Processing, vol. 59, no. 4, pp. 1409-1420, Apr. 2011.
%[2] M. Feldmann, "Tracking von Objektgruppen und ausgedehnten
%    Zielobjekten," Ph.D. dissertation, KIT-Fakultät für Informatik des
%    Karlsruher Instituts für Technologie, Karlsruhe, Germany, 30 Nov.
%    2018.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function vals=PDF(wVals,nu,lambda)
%%PDF Evaluate the PDF of the distribution at one or more points.
%
%INPUTS: wVals A dXN matrix of the sets of eigenvalues where it is desired
%              that the PDF be evaluated. wVals(:,k) is the kth set of
%              (unordered) eigenvalues. It is required that d>1.
%           nu The number of degrees of freedom of the distribution. Note
%              that nu>=d.
%       lambda The scale matrix of the underlying Wishart distribution is
%              assumed to be lambda*eye(d,d). lambda>0. If omitted or an
%              empty matrix is passed, the default is 1.
%
%OUTPUTS: vals The NX1 set of PDF values evaluated at the points in wVals.
%
%The PDF is given in Corollary 3.2.19 in [1], which is a simplification of
%a more general formula, originally from Equation 58 in [2]. The formula in
%[2] has a normalization constant that is difficult to evaluate. The
%formula of [1] is implemented in this function.
%
%A simplification for the 2D case (d=2) arises in [2], [3] when considering
%displaying the extension of the target when performing extended target
%tracking. The expressions in [3] are corrected compared to those in [2].
%The simplification of [3] is implemented for the d=2 case to be more
%efficient.
%
%REFERENCES:
%[1] R. J. Muirhead, Aspects of Multivariate Statistical Theory. Hoboken:
%    John Wiley & Sons, 2005.
%[2] A. T. James, "Distributions of matrix variates and latent roots
%    derived from normal samples," Annals of Mathematical Statistics, vol.
%    35, no. 2, pp. 475-501, Jun. 1964.
%[3] M. Feldmann, D. Fränken, and W. Koch, "Tracking of extended objects
%    and group targets using random matrices," IEEE Transactions on Signal
%    Processing, vol. 59, no. 4, pp. 1409-1420, Apr. 2011.
%[4] M. Feldmann, "Tracking von Objektgruppen und ausgedehnten
%    Zielobjekten," Ph.D. dissertation, KIT-Fakultät für Informatik des
%    Karlsruher Instituts für Technologie, Karlsruhe, Germany, 30 Nov.
%    2018.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.   

if(nargin<3||isempty(lambda))
    lambda=1;
end

d=size(wVals,1);
N=size(wVals,2);

if(d<2)
    error('The distribution must be at least 2D.')
end
vals=zeros(N,1);

if(d==2)
    for k=1:N
        w=wVals(:,k);
        w=sort(w,'descend');
        %The 2D distribution.
        w1=w(1);
        w2=w(2);

        %Below is equivalent to:
        %val=exp(-(1/(2*lambda))*(w1+w2))*(w1-w2)*(w1*w2)^((n-3)/2)/(lambda^n*4*gamma(n-1));
        vals(k)=exp(-(1/(2*lambda))*(w1+w2)+log(w1-w2)+((nu-3)/2)*(log(w1)+log(w2))-nu*log(lambda)-log(4)-gammaln(nu-1));
    end
else
    for k=1:N
        w=wVals(:,k);
        w=sort(w,'descend');
        term1=(d^2/2)*log(pi)-((d*nu/2)*log(2*lambda)+gammalnMultiDim(d/2,d)+gammalnMultiDim(nu/2,d))-(1/(2*lambda))*sum(w)+((nu-d-1)/2)*sum(log(w));

        term2=0;
        for j=2:d
            for i=1:(j-1)
                term2=term2+log(w(i)-w(j));
            end
        end

        vals(k)=exp(term1+term2);
    end
end
end

function val=maxMinEigRegion2D(lMin,lMax,nu,lambda)
%%MAXMINEIGREGION2D Determine Pr{lMin<=eigMin,eigMax<=lMax), where eigMin
%            and eigMax are the minimum and maximum eigenvalues of a sample
%            of the 2D Wishart distribution with nu degrees of freedom and
%            scale matrix lambda*eye(2,2).
%
%INPUTS: lMin, lMax The scalar lower bound on the minimum eigenvalue and
%          the scalar upper bound on the maximum eigenvalue. lMin>=0 and
%          lMax<=Inf.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%   lambda The scale matrix of the underlying Wishart distribution is
%          assumed to be lambda*eye(d,d). lambda>0. If omitted or an empty
%          matrix is passed, the default is 1.
%
%OUTPUTS: vals A matrix of probabiltiies having a dimensionality consistent
%              with lMinVals and lMaxVals.
%
%The 2D PDF is given in Appendix D.2 of [1]. The difficult integral is
%solved using the expression in terms of incomplete gamma functions that is
%given in Appendix D.1 of [1]. The expressions omit lambda, but the
%incorporation of lambda just scales lMin and lMax in the computation.
%
%EXAMPLE:
%Here, we demonstrate that the probability in the region between the
%bounded minimum and maximum eigenvalues is consistent with one found via
%Monte Carlo simulation.
% lambda=21.2;
% A=lambda*eye(2,2);
% nu=6;
% numMC=2e5;
% lMin=100;
% lMax=200;
% numInRegion=0;
% for k=1:numMC
%     eigVals=eig(WishartD.rand(A,nu));
%     maxEig=max(eigVals);
%     minEig=min(eigVals);
%     
%     if(lMin<=minEig&&maxEig<=lMax)
%         numInRegion=numInRegion+1;
%     end
% end
% empVal=numInRegion/numMC
% exactVal=StdWishartEigenvalueD.maxMinEigRegion2D(lMin,lMax,nu,lambda)
%Both solutions will tend to be around 0.556, indicating consistency of the
%techniques.
%
%REFERENCES:
%[1] M. Feldmann, "Tracking von Objektgruppen und ausgedehnten
%    Zielobjekten," Ph.D. dissertation, KIT-Fakultät für Informatik des
%    Karlsruher Instituts für Technologie, Karlsruhe, Germany, 30 Nov.
%    2018.
%[2] R. J. Muirhead, Aspects of Multivariate Statistical Theory. Hoboken:
%    John Wiley & Sons, 2005.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(lMin>lMax)
        error('lMin must be <=lMax.')
    end

    if(nargin<3||isempty(lambda))
        lambda=1;
    end

    lMin=lMin/lambda;
    lMax=lMax/lambda;
    
    beta=1+(1/2)*(nu-3);
    coeff=-sqrt(pi)*exp(((1/2)-(nu/2))*log(2)-gammaln(nu/2));
    
    expTerm1=exp(-lMax/2);
    expTerm2=lMax.^beta;
    
    if(~isfinite(lMax)||~isfinite(expTerm1*expTerm2))
        %Use the limit for infinite lMax (or lMax large enough to cause a
        %0*Inf.
        val=1-gammainc(lMin,nu-1)+...
            coeff*exp(-lMin/2)*lMin^beta.*(1-gammainc(lMin/2,beta));
    else
        val=gammainc(lMax,nu-1)-gammainc(lMin,nu-1)+...
            coeff*(expTerm1.*expTerm2+exp(-lMin/2).*lMin.^beta).*...
            (gammainc(lMax/2,beta)-gammainc(lMin/2,beta));
    end

    %Deal with possible finite precision issues.
    val=max(0,val);
    val=min(1,val);
end


function vals=CDF2DMax(x,nu,lambda)
%%CDF2DMAX This is the CDF of the marginal PDF of the maximum eigenvalue of
%          a 2D random variable generated by a Wishart distribution with nu
%          degrees of freedom and scale matrix lambda*eye(2,2).
%
%INPUTS: x A matrix of points at which CDF values are desired.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=d.
%   lambda The scale matrix of the underlying Wishart distribution is
%          assumed to be lambda*eye(d,d). lambda>0. If omitted or an empty
%          matrix is passed, the default is 1.
%
%OUTPUTS: vals A matrix of CDF values have the same dimensions as x.
%
%The CDF is explicitely given in Appendix D.2 of [1]. It is based on
%integrating the PDF given in Corollary 3.2.19 in [2].
%
%EXAMPLE:
%Here, we demonstrate that the CDF of the maximum eigenvalue of the Wishart
%distribution is consistent with one found via Monte Carlo simulation and
%the formation of an empirical CDF value.
% lambda=21.2;
% A=lambda*eye(2,2);
% nu=6;
% numMC=2e5;
% x=200;
% numInRegion=0;
% for k=1:numMC
%     maxEig=max(eig(WishartD.rand(A,nu)));
%     if(maxEig<=x)
%         numInRegion=numInRegion+1;
%     end
% end
% empVal=numInRegion/numMC
% exactVal=StdWishartEigenvalueD.CDF2DMax(x,nu,lambda)
%Both solutions will tend to be around 0.61, indicating consistency of the
%techniques.
%
%REFERENCES:
%[1] M. Feldmann, "Tracking von Objektgruppen und ausgedehnten
%    Zielobjekten," Ph.D. dissertation, KIT-Fakultät für Informatik des
%    Karlsruher Instituts für Technologie, Karlsruhe, Germany, 30 Nov.
%    2018.
%[2] R. J. Muirhead, Aspects of Multivariate Statistical Theory. Hoboken:
%    John Wiley & Sons, 2005.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(lambda))
        lambda=1;
    end

    x=x/lambda;

    vals=gammainc(x,nu)-exp(-x/2-gammaln(nu/2)+(1/2)*(nu-1)*log(x/2)+(1/2)*log(pi)).*gammainc(x/2,(nu+1)/2);

    %Deal with possible finite precision issues.
    vals=max(0,vals);
    vals=min(1,vals);
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
