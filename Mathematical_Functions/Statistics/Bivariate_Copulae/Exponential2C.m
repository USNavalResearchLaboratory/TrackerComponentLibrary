classdef Exponential2C
%%EXPONENTIAL2C A collection of functions for the bivariate Exponential
%               copulae. Implemented functions include: PDF, CDF.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%[2] X. Zeng, J. Ren, Z. Wang, S. Marshall, and T. Durrani, "Copulas for
%    statistical signal processing (part i): Extensions and
%    generalization," Signal Processing, vol. 94, pp. 691-702, 2014,
%    ISSN: 0165-1684. DOI: https://doi.org/10.1016/j.sigpro.2013.07.009.
%    [Online]. Available: http://www.sciencedirect.com/science/article/pii/S0165168413002880.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function jprob = PDF(u,rho)
            %PDF Generates a joint probability for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % rho: A scalar coupling parameter in the interval
            %      (0,1). Defaults to 0.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %EXAMPLE 1: Make a contour plot of the density.
            %x = linspace(0,1,1e2);
            %[X,Y] = meshgrid(x);
            %Z = zeros(length(X),length(Y));
            %for idx = 1:length(X)
            %    for iidx = 1:length(Y)
            %        Z(idx,iidx) = Exponential2C.PDF([X(idx,iidx);Y(idx,iidx)],0.5);
            %    end
            %end
            %contourf(X,Y,Z)
            %
            %EXAMPLE 2: Generate an exact joint density plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            % for iidx = 1:length(Y)
            % Z(idx,iidx) = Exponential2C.PDF([X(idx,iidx);Y(idx,iidx)],0.5);
            % X(idx,iidx) = GumbelD.invCDF(X(idx,iidx),1,2);
            % Y(idx,iidx) = GaussianD.invCDF(Y(idx,iidx));
            % Z(idx,iidx) = Z(idx,iidx).*GumbelD.PDF(X(idx,iidx),1,2).*GaussianD.PDF(Y(idx,iidx));
            % end
            % end
            % surf(X,Y,Z)
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('u','var') || isempty(u)
                jprob = [];
                return
            end
            if ~exist('rho','var') || isempty(rho)
                rho = 0;
            elseif rho<0 || rho>=1
                error('Correlation coefficient must be in the interval [0,1).')
            end
            
            a = log(1-u(1,:));
            b = log(1-u(2,:));
            c = 2*sqrt(rho*a*b)/(1-rho);
            
            jprob = exp(rho*(a+b)/(1-rho)).*besseli(0,c)/(1-rho);
            
            %NaN results from 0*inf in the calculation.
            %These should be 0s.
            nanIdxs = isnan(jprob);
            jprob(nanIdxs) = 0;
        end
        
        function jdist = CDF(u,rho,lambdas)
            %CDF Generates a joint distribution for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % rho: A scalar coupling parameter in the interval
            %      (0,1). Defaults to 0.
            % lambdas: A 2-by-1 or 1-by-2 vector containing the two rate
            %          constants lambda1 and lambda2 associated with
            %          dimension 1 and dimension 2 respectively.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint distribution values.
            %
            %EXAMPLE 1: Make a contour plot of the distribution.
            %x = linspace(0,1,1e2);
            %[X,Y] = meshgrid(x);
            %Z = zeros(length(X),length(Y));
            %for idx = 1:length(X)
            %    for iidx = 1:length(Y)
            %        Z(idx,iidx) = Exponential2C.CDF([X(idx,iidx);Y(idx,iidx)],0.5);
            %    end
            %end
            %contourf(X,Y,Z)
            %
            %EXAMPLE 2: Generate an exact joint distribution plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            % for iidx = 1:length(Y)
            % Z(idx,iidx) = Exponential2C.CDF([X(idx,iidx);Y(idx,iidx)],0.5);
            % X(idx,iidx) = GumbelD.invCDF(X(idx,iidx),1,2);
            % Y(idx,iidx) = GaussianD.invCDF(Y(idx,iidx));
            % Z(idx,iidx) = Z(idx,iidx);
            % end
            % end
            % surf(X,Y,Z)
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('u','var') || isempty(u)
                jdist = [];
                return
            end
            if ~exist('rho','var') || isempty(rho)
                rho = 0;
            elseif rho<0 || rho>=1
                error('Correlation coefficient must be in the interval [0,1).')
            end
            if ~exist('lambdas','var') || isempty(lambdas)
                lambda1 = 1;
                lambda2 = 1;
            else
                lambda1 = lambdas(1);
                lambda2 = lambdas(2);
            end
            
            a = -log(1-u(1,:))/lambda1;
            b = -log(1-u(2,:))/lambda2;
            c = 2*sqrt(rho*lambda1*lambda2*a*b)/(1-rho);
            
            f = @(x,y) (lambda1*lambda2/(1-rho))*exp(-(lambda1*a+lambda2*b)/(1-rho))...
                       *besselj(0,c);
            
            for idx = 1:length(a)
                for iidx = 1:length(b)
                    jdist = integral2(f,0,alim,0,blim);
                end
            end
            
            %The CDF takes the max of (a+b-1).^(-1/theta) and 0.
            negIdxs = jdist<0;
            jdist(negIdxs) = 0;
            
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
