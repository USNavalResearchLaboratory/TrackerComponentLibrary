classdef PearsonIID
%%PEARSONIID Functions to handle the multivariate Pearson type II
%            distribution.
%Implemented methods are: mean, PDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static) 
function val=mean(d)
%%MEAN Obtain the mean of the multivariate Pearson-II distribution.
%
%INPUTS: d The dimensionality of the distribution.
%
%OUTPUTS: val The dX1 mean of the distribution. This is just zero, because,
%             as can be seen in [1], the distribution is symmetric about
%             the origin.
%   
%REFERENCES:
%[1] M. Khalafi and M. Azimmohseni, "Multivariate Pearson type II
%    distribution: Statistical and mathematical features," Probability and
%    Mathematical Statistics, vol. 34, no. 1, pp. 119-126, 2014.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %The distribution is symmetric about the origin.
    val=zeros(d,1);
end
    
function val=PDF(x,a)
%%PDF Evaluate the PDF of a multivariate Pearson type II distribution. The
%     distribution is given in [1] and, as can be seen in Chapter 4.2.1 of
%     [2], when a=2, the distribution is the same as the Epanechnikov
%     kernel.
%
%INPUTS: x The dXN set of N d-dimensional points at which values of the
%          distribution are desired.
%        a The shape parameter of the distribution. If a=2, then one has
%          the Epanechnikov kernel. The default if omitted or an empty
%          matrix is passed is 2.
%
%OUTPUTS: val The 1XN set of values of the distribution at the points in x.
%   
%REFERENCES:
%[1] M. Khalafi and M. Azimmohseni, "Multivariate Pearson type II
%    distribution: Statistical and mathematical features," Probability and
%    Mathematical Statistics, vol. 34, no. 1, pp. 119-126, 2014.
%[2] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(a))
        a=2;%Epanechnikov kernel
    end
    
    d=size(x,1);%The dimensionality of the points.
    Ns=size(x,2);
    val=zeros(1,Ns);

    const=exp(gammaln(a+d/2)-gammaln(a))*pi^(-d/2);

    xMags2=sum(x.*x,1);
    sel=xMags2<1;
    val(sel)=const*(1-xMags2(sel)).^(a-1);
end
    
function vals=rand(N,d,a)
%%RAND Generate random samples of the multivariate Pearson II distribution.
%
%INPUTS: N the scalar number of samples to generate.
%        d The dimensionality of the samples.
%        a The shape parameter of the distribution. If a=2, then one has
%          the Epanechnikov kernel. The default if omitted or an empty
%          matrix is passed is 2.
%
%OUTPUTS: vals A dXN set of N random samples of the distribution.
%
%The method for sampling using the uniform and beta distributions is given
%in Theorem 2.1 in [1].
%
%EXAMPLE:
%Here, we show that in 1D, the histrogram of the points matches the PDF
%values.
% N=1e4;
% numPlotPoints=500;
% d=1;
% a=50;
% 
% ySamp=PearsonIID.rand(N,d,a);
% 
% x=linspace(-1,1,numPlotPoints);
% y=PearsonIID.PDF(x,a);
% 
% figure(1)
% clf
% hold on
% histogram(ySamp,'Normalization','pdf')
% plot(x,y,'linewidth',2);
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] M. Khalafi and M. Azimmohseni, "Multivariate Pearson type II
%    distribution: Statistical and mathematical features," Probability and
%    Mathematical Statistics, vol. 34, no. 1, pp. 119-126, 2014.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(a))
        a=2;%Epanechnikov kernel
    end
    
    V=BetaD.rand([1,N],d/2,a);
    %d-dimensional random vectors uniformly distributed on the unit sphere.
    U=randDirVec(d,N);
    
    vals=bsxfun(@times,U,sqrt(V));
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
