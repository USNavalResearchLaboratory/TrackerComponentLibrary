classdef InverseWishartD
%%INVERSEWISHARTD Functions to handle the inverse Wishart distribution.
%Implemented methods are: PDF, mean,rand
%
%DEPENDENCIES: WishartD.m
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)
    
function val=mean(Psi,nu)
%%MEAN Find the mean of a particular inverse Wishart distribution
%
%INPUTS:   Psi  The D X D positive-definite, symmetric precision matrix of
%               the distribution.
%           nu  The number of degrees of freedom of the distribution. Note
%               that nu>=D.
%
%OUTPUTS: val  The DXD mean of the inverse Wishart distribution.
%
%The inverse Wishard distribution is discussed in Chapter 3 of [1].
%
%REFERENCES:
%[1] K. V. Mardia, J. T. Kent, and J. M. Bibby, Multivariate Analysis.
%    Academic Press., 1979
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

D=size(Psi,1);
val=Psi/(nu-D-1);

end
    
function val=PDF(X,Psi,nu)
%%PDF Evaluate the inverse Wishart distribution for a particular
%     matrix.
%
%INPUTS:    X   The D X D symmetric matrix at which the PDF of the
%               distribution is to be evaluated.
%          Psi  The D X D positive-definite, symmetric precision matrix of
%               the distribution.
%           nu  The number of degrees of freedom of the distribution. Note
%               that nu>=D.
%
%OUTPUTS:   val The value of the PDF of the inverse Wishart distribution at
%               X with the given parameters.
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
    
function val=rand(Psi,nu)
%%RAND      Generate an inverse  Wishart random matrix with the given number of
%           degrees of freedom and precision matrix.
%
%INPUTS:   Psi  The D X D positive-definite, symmetric precision matrix of
%               the distribution.
%           nu  The number of degrees of freedom of the distribution. Note
%               that nu>=D.
%
%OUTPUTS: X  A DXD inverse Wishart random matrix.
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
    
    val=WishartD.rand(inv(Psi),nu);
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
