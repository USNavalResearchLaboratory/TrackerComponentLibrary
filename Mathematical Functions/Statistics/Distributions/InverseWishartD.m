classdef InverseWishartD
%%INVERSEWISHARTD Functions to handle the inverse Wishart distribution
%                 (also known as the inverted Wishart distribution).
%                 Inverse Wishart distributions are used in some extended
%                 target tracking algorithms that model the target extent
%                 as an inverse Wishart distributed random variable.
%Implemented methods are: PDF, mean, var, cov, secondMatMoment,
%                         secondMoments, rand
%
%DEPENDENCIES: WishartD.m
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)
function val=mean(Psi,nu)
%%MEAN Find the mean of a particular inverse Wishart distribution
%
%INPUTS: Psi The dXd positive-definite, symmetric precision matrix of
%            the distribution.
%         nu The number of degrees of freedom of the distribution. Note
%            hat nu>d+1. 
%
%OUTPUTS: val The dXd mean of the inverse Wishart distribution.
%
%The inverse Wishard distribution is discussed in Chapter 3 of [1].
%
%EXAMPLE:
%Here, we show that the mean is consistent with a mean obtained via Monte
%Carlo simulation.
% Psi=inv([55,  7, 12, 17;
%           7, 45, 17, 22;
%          12, 17, 35, 27;
%          17, 22, 27, 25]);
% nu=22;
% numMC=1e5;
% XMean=zeros(4,4);
% XVals=zeros(4,4,numMC);
% for k=1:numMC
%     X=InverseWishartD.rand(Psi,nu);
%     XMean=XMean+X;
%     XVals(:,:,k)=X;
% end
% XMean=XMean/numMC;
% XMeanExact=InverseWishartD.mean(Psi,nu);
% RelErr=max(abs((XMean(:)-XMeanExact(:))./XMeanExact(:)))
%The relative error tends to be around 0.002, indicating good agreement.
%
%REFERENCES:
%[1] K. V. Mardia, J. T. Kent, and J. M. Bibby, Multivariate Analysis.
%    Academic Press., 1979
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

d=size(Psi,1);
val=Psi/(nu-d-1);
end
    
function val=PDF(X,Psi,nu)
%%PDF Evaluate the inverse Wishart distribution for a particular
%     matrix.
%
%INPUTS: X The D X D symmetric matrix at which the PDF of the distribution
%          is to be evaluated.
%      Psi The D X D positive-definite, symmetric precision matrix of the
%          distribution.
%       nu The number of degrees of freedom of the distribution. Note that
%          nu>=D.
%
%OUTPUTS: val The value of the PDF of the inverse Wishart distribution at
%             X with the given parameters.
%
%The inverse Wishart distribution is discussed in Chapter 7.7 of [1].
%
%REFERENCES:
%[1] T. W. Anderson, An Introduction to Multivariate Statistical Analysis,
%    3rd ed. Hoboken, NJ: Wiley-Interscience, 2003.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    D=size(X,1);

    if(rank(X)<D)
        val=0;
        return;
    end
    
    %The normalization constant
    B=det(Psi)^(nu/2)/(2^(nu*D/2)*pi^(D*(D-1)/4)*prod(gamma((nu+1-(1:D))/2)));
    val=B*det(X)^(-(nu+D+1)/2)*exp(-0.5*trace(Psi/X));
end

function val=var(Psi,nu)
%%VAR Obtain the variance of the inverse Wishart distribution. Though the
%     random variable is a matrix, the variance is E{(X-xMean)*(X-XMean)},
%     where xMean is the mean of the distribution.
%
%INPUTS: Psi The dXd positive-definite, symmetric precision matrix of
%            the distribution.
%         nu The number of degrees of freedom of the distribution. Note
%            that nu>d+3 for a finite solution to exist.
%
%OUTPUS: val The dXd variance matrix.
%
%The variance is just computed from the identity that the variance is
%E{X^2}-E{X}^2.
%
%EXAMPLE:
%Here, we show that the variance is consistent with a variance obtained via
%Monte Carlo simulation.
% Psi=inv([55,  7, 12, 17;
%           7, 45, 17, 22;
%          12, 17, 35, 27;
%          17, 22, 27, 25]);
% nu=40;
% numMC=1e5;
% XMean=zeros(4,4);
% XVals=zeros(4,4,numMC);
% for k=1:numMC
%     X=InverseWishartD.rand(Psi,nu);
%     XMean=XMean+X;
%     XVals(:,:,k)=X;
% end
% XMean=XMean/numMC;
% 
% XVar=zeros(4,4);
% for k=1:numMC 
%     diff=XVals(:,:,k)-XMean;
%     XVar=XVar+diff*diff;
% end
% XVar=XVar/numMC;%Sample variance
% XVarTrue=InverseWishartD.var(Psi,nu);
% RelErr=max(abs((XVar(:)-XVarTrue(:))./XVarTrue(:)))
%The relative error tends to be around than 0.003, which is good.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    EX2=InverseWishartD.secondMatMoment(Psi,nu);
    EX=InverseWishartD.mean(Psi,nu);
    val=EX2-EX*EX;
end

function matVal=cov(Psi,nu)
%%COV Calculate the covariance of all of the elements of the dXd inverse
%     Wishart matrix with each other. Since the inverse Wishart
%     distribution is a matrix distribution, the covariance is stored in a
%     tensor (4D matrix).
%
%INPUTS: Psi The dXd positive-definite, symmetric precision matrix of
%            the distribution.
%         nu The number of degrees of freedom of the distribution. Note
%            that nu>d+3 for a finite solution to exist.
%
%OUTPUTS: matVal The dXdXdXd matrix of covariance values
%                matVal(i1,i2,i3,i4) is the covariance of X(i1,i2) with
%                X(i3,i4).
%
%An expression for the covariance, though with a different
%parameterization, is given on page 113 of [1] and has been adapted for use
%here.
%
%REFERENCES:
%[1] A. K. Gupta and D. K. Nagar, Matrix Variate Distributions, London:
%    Chapman & Hall/CRC, 1999.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(Psi,1);
    matVal=zeros(d,d,d,d);
    
    for i=1:d
        for j=i:d 
            for k=1:d
                for l=k:d
                    val=((nu-d-1)*(Psi(i,k)*Psi(j,l)+Psi(i,l)*Psi(k,j))+2*Psi(i,j)*Psi(k,l))/((nu-d)*(nu-d-1)^2*(nu-d-3));
                    matVal(i,j,k,l)=val;
                    matVal(j,i,k,l)=val;
                    matVal(i,j,l,k)=val;
                    matVal(j,i,l,k)=val;
                end
            end
        end
    end
end


function matVal=secondMatMoment(Psi,nu)
%%SECONDMATMOMENT Compute the second matrix moment of the inverse Wishart
%                 distribution. That is E{X*X}.
%
%INPUTS: Psi The dXd positive-definite, symmetric precision matrix of
%            the distribution.
%         nu The number of degrees of freedom of the distribution. Note
%            that nu>d+3 for a finite solution to exist.
%
%OUTPUS: matVal The dXd second moment matrix.
%
%The second moment is given in [1].
%
%EXAMPLE:
%Here, we show that the second moment is consistent with a second moment
%obtained via Monte Carlo simulation.
% Psi=inv([55,  7, 12, 17;
%           7, 45, 17, 22;
%          12, 17, 35, 27;
%          17, 22, 27, 25]);
% nu=40;
% numMC=1e5;
% EX2=zeros(4,4);
% for k=1:numMC
%     X=InverseWishartD.rand(Psi,nu);
%     EX2=EX2+X*X;
% end
% EX2=EX2/numMC;
% EX2Exact=InverseWishartD.secondMatMoment(Psi,nu);
% RelErr=max(abs((EX2(:)-EX2Exact(:))./EX2Exact(:)))
%The relative error tends to be around 0.003, which is good.
%
%REFERENCES:
%[1] L. R. Haff, "Identities for the inverse Wishart distribution with
%    computational results in linear and quadratic discrimination,"
%    Sankhya: The Indian Journal of Statistics, Series B, vol. 44, no. 3,
%    pp. 245-258, Dec. 1982.
%
%April 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(Psi,1);
    matVal=(trace(Psi)*Psi+(nu-d-1)*Psi*Psi)/((nu-d)*(nu-d-1)*(nu-d-3));
end

function matVal=secondMoments(Psi,nu)
%%SECONDMOMENTS Calculate the second moment (and all cross moments) of all
%     of the elements of the dXd inverse Wishart distribution. Since the
%     inverse Wishart distribution is a matrix distribution, the second
%     moments are stored in a tensor (4D matrix).
%   
%INPUTS: Psi The dXd positive-definite, symmetric precision matrix of
%            the distribution.
%         nu The number of degrees of freedom of the distribution. Note
%            that nu>d+3 for a finite solution to exist.
%
%OUTPUTS: matVal The dXdXdXd matrix of second moments matVal(i1,i2,i3,i4)
%                is the expected value E(X(i1,i2)*X(i3,i4)).
%
%The second moments are found using the identity E{X*Y}=Cov{X,Y}+E{X}*E{Y}
%and the cov and mean methods of this class.
%
%EXAMPLE:
%Here, we compare the second-order moments obtained from Monte Carlo runs
%to those obtained from this function.
% Psi=[55,  7, 12, 17;
%       7, 45, 17, 22;
%      12, 17, 35, 27;
%      17, 22, 27, 25];
% nu=40;
% d=size(Psi,1);
% 
% m2Exact=InverseWishartD.secondMoments(Psi,nu);
% 
% numMonteCarlo=1e5;
% m2Mat=zeros(4,4,4,4);
% for curRun=1:numMonteCarlo
%     X=InverseWishartD.rand(Psi,nu);
%     
%     %Find all second moments.
%     for i=1:d
%         for j=1:d 
%             for k=1:d
%                 for l=1:d
%                     m2Mat(i,j,k,l)=m2Mat(i,j,k,l)+X(i,j)*X(k,l);
%                 end
%             end
%         end
%     end
% end
% m2Mat=m2Mat/numMonteCarlo;
% %Compare to the analytic value.
% max(abs(m2Mat(:)-m2Exact(:))./abs(m2Exact(:)))
%The values should agree within about half a percent. Note that decreasing
%nu increases the number fo Monte Carlo runs needed for an accurate result.   
%   
%May 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    d=size(Psi,1);

    meanVal=InverseWishartD.mean(Psi,nu);
    
    prodVals=zeros(d,d,d,d);
    for i=1:d
        for j=1:d
            for k=1:d
                for l=1:d
                    prodVals(i,j,k,l)=meanVal(i,j)*meanVal(k,l);
                end
            end
        end
    end
    covMat=InverseWishartD.cov(Psi,nu);
 
    matVal=covMat+prodVals;
end

function val=rand(Psi,nu)
%%RAND Generate an inverse  Wishart random matrix with the given number of
%      degrees of freedom and precision matrix.
%
%INPUTS: Psi The dXd positive-definite, symmetric precision matrix of the
%            distribution.
%         nu The number of degrees of freedom of the distribution. Note
%            that nu>=d. This must be an integer.
%
%OUTPUTS: X A dXd inverse Wishart random matrix.
%
%The inverse Wishart distribution is discussed in Chapter 7.7 of [1]. It is
%related to a transformed Wishard random matrix. The transformation is used
%here.
%
%REFERENCES:
%[1] T. W. Anderson, An Introduction to Multivariate Statistical Analysis,
%    3rd ed. Hoboken, NJ: Wiley-Interscience, 2003.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    val=pinv(WishartD.rand(inv(Psi),nu));
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
