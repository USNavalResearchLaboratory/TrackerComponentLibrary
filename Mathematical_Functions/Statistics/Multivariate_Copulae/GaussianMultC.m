classdef GaussianMultC
%%GAUSSIANMULTC A collection of functions for the Gaussian multivariate
%               copulae. Implemented functions are: PDF.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function jprob = PDF(u,R)
        %PDF Generates a joint probability for column vector u or a row
        %    vector of joint probabilities for matrix u with each column as
        %    the input vector.
        %
        %INPUTS: u: A d-by-n matrix where each column represents a
        %           d-dimensional random vector and n is the number of
        %           random vectors. Entries must be in [0,1].
        %        R: A d-by-d correlation matrix.
        %
        %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        if ~exist('u','var') || isempty(u)
            jprob = [];
            return
        end

        s2 = size(u,2);
        I = eye(size(R));
        jprob = zeros(1,s2);

        for col = 1:s2
            v = GaussianD.invCDF(u(:,col));
            jprob(col) = exp(-0.5*v'*(inv(R)-I)*v)/sqrt(det(R));
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
