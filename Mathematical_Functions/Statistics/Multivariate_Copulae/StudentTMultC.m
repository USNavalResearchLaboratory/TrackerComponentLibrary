classdef StudentTMultC
%%STUDENTMULTC A collection of functions for the Student t-copulae.
%              Implemented functions include: PDF.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function jprob = PDF(u,sigma,nu)
            %PDF Generates a joint probability for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: u: A d-by-n matrix where each column represents a
            %           d-dimensional random vector and n is the number of
            %           random vectors. Entries must be in [0,1].
            %    sigma: A d-by-d matrix defining the scaling of the
            %           distribution. Sigma is proportional to the
            %           covariance with a proportionality constant of
            %           nu/(nu-2) for nu>2.
            %       nu: Scalar degrees of freedom.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('u','var') || isempty(u)
                jprob = [];
                return
            end
            
            [s1,s2] = size(u);
            
            if ~exist('sigma','var') || isempty(sigma)
                sigma = eye(s1,s1);
            end
            
            rho = sigma;
            for idx = 1:s1
                for iidx = 1:s1
                    rho(idx,iidx) = rho(idx,iidx)/sqrt(sigma(idx)*sigma(iidx));
                end
            end
            jprob = ones(1,s2)*(1/sqrt(det(rho)))*(gamma((nu+s1)/2)/gamma(nu/2))*(gamma(nu/2)/gamma((nu+1)/2))^s1;

            for col = 1:s2
                v = u(:,col);
                jprob(col) = jprob(col)*prod((1+StudentTD.invCDF(v).^2/nu).^((nu+1)/2))/(1+GaussianD.invCDF(v)'*pinv(rho)*GaussianD.invCDF(v)/nu)^((nu+s1)/2);
            end
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
