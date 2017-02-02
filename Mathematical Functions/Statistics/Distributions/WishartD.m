classdef WishartD
%Functions to handle the Wishart distribution.
%Implemented methods are: PDF, mean, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(A,nu)
%%MEAN Find the mean of a particular Wishart distribution
%
%INPUTS:    A   The D X D positive-definite, symmetric scale matrix of the
%               distribution.
%           nu  The number of degrees of freedom of the distribution. Note
%               that nu>=D.
%
%OUTPUTS: val  The DXD mean of the Wishart distribution.
%
%The Wishart distribution and its generation are discussed in chapters 5
%and 8 of [1].
%
%REFERENCES:
%[1] M. L. Eaton, Multivariate Statistics: A Vector Space Approach, ser.
%    Institute of Mathematical Statistics Lecture Notes-Monograph Series.
%    Beachwood, Ohio: Institute of Mathematical Statistics, 2007, vol. 53.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=A*nu;

end
    
function val=PDF(X,A,nu)
%%PDF Evaluate the Wishart distribution for a particular matrix.
%
%INPUTS:    X   The D X D symmetric matrix at which the PDF of the
%               distribution is to be evaluated.
%           A   The D X D positive-definite, symmetric scale matrix of the
%               distribution.
%           nu  The number of degrees of freedom of the distribution. Note
%               that nu>=D.
%
%OUTPUTS:   val The value of the PDF of the Wishart distribution at X with
%               the given parameters.
%
%An expression for this Wishart distribution with the normalization
%constant is given in Chapter 2.3.6 of [1].
%
%REFERENCES:
%[1] C. M. Bishop, Pattern Recognition and Machine Learning. Cambridge,
%    United Kingdom: Springer, 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    D=size(X,1);
    
    if(rank(X)<D)
        val=0;
        return;
    end
    
    %The normalization constant
    B=det(A)^(-nu/2)/(2^(nu*D/2)*pi^(D*(D-1)/4)*prod(gamma((nu+1-(1:D))/2)));

    val=B*det(X)^((nu-D-1)/2)*exp(-0.5*trace(A\X));
end
    
function X=rand(A,nu)
%%RAND      Generate a Wishart random matrix with the given number of
%           degrees of freedom and scale matrix.
%
%INPUTS:    A   The D X D positive-definite, symmetric scale matrix of the
%               distribution.
%           nu  The number of degrees of freedom of the distribution. Note
%               that nu>=D.
%
%OUTPUTS: X  A DXD Wishart random matrix.
%
%The Wishart distribution and its generation are discussed in chapters 5
%and 8 of [1]. The Wishart distribution is just the outer product of DXnu
%normally generated random matrices, where each row is generated having
%zero-mean and covariance matrix A.
%
%REFERENCES:
%[1] M. L. Eaton, Multivariate Statistics: A Vector Space Approach, ser.
%    Institute of Mathematical Statistics Lecture Notes-Monograph Series.
%    Beachwood, Ohio: Institute of Mathematical Statistics, 2007, vol. 53.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    D=size(A,1);

    S=chol(A,'lower')*randn(D,nu);
    X=S*S';
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
