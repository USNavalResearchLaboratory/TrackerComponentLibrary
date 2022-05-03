classdef Empirical2C
%%EMPIRICAL2C A collection of functions for the bivariate empirical copula.
%             Implemented functions include: PDF and CDF.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function [jprob,points] = PDF(u,N)
            %PDF Generates a joint copula density approximation on a
            %    grid of the unit square given marginal CDF values for
            %    a collection of realized values of a bivariate random
            %    vector.
            %
            %INPUT:
            % u: A 2-by-numpoints matrix where each column represents a
            %    2-dimensional vector of the marginal CDFs for each random
            %    variable's realization and numpoints is the number
            %    of random vectors. Entries must be in [0,1].
            % N: The number of grid points to use for dividing each axis
            %    of the unit square. The resulting grid of points will be
            %    of size N-by-N. Defaults to 100.
            %
            %OUTPUT:
            % jprob: A 1-by-N^2 vector of joint probabilities.
            % points: A 2-by-N^2 vector of points corresponding to the
            %         probabilities in jprob.
            %
            %EXAMPLE 1: Generates a linearly correlated vector of Gaussian
            %           marginals and plots the copula density.
            % % Generate some data
            % numd1 = 1e5;
            % data = randn(2,numd1);
            % sigma = chol([1,0.5;0.5,1],'lower');
            % data = sigma*data;
            % d1 = data(1,:);
            % d2 = data(2,:);
            % u1 = GaussianD.CDF(d1);
            % u2 = GaussianD.CDF(d2);
            % [jprob,points] = Empirical2C.PDF([u1;u2],50);
            % X = reshape(points(1,:),50,50);
            % Y = reshape(points(2,:),50,50);
            % surf(X,Y,reshape(jprob,50,50))
            %
            %December 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if(nargin<2||isempty(N))
                % Generate a uniform grid of (0,1)x(0,1), 100x100 points
                N=100;
            end

            x = linspace(0,1,N);
            [X,Y] = meshgrid(x);
            h1 = x(2)-x(1);
            h2 = h1;

            u1 = u(1,:);
            u2 = u(2,:);

            % Compute copula distribution values
            points = [X(:),Y(:)]';
            numpoints = size(points,2);
            counts = zeros(1,numpoints);
            for idx = 1:numpoints
                counts(idx) = sum((u1<points(1,idx))&(u2<points(2,idx)));
            end
            counts = counts/sum(counts);

            % Compute differences to get a density
            jprob = reshape(counts,size(X,2),size(Y,1));
            jprob = gradient(jprob)./h1;
            jprob = gradient(jprob')'./h2;
            jprob = jprob(:);
        end
        
        function [jdist,points] = CDF(u,N)
            %CDF Generates a joint copula distribution approximation on a
            %    grid of the unit square given marginal CDF values for
            %    a collection of realized values of a bivariate random
            %    vector.
            %
            %INPUT:
            % u: A 2-by-numpoints matrix where each column represents a
            %    2-dimensional vector of the marginal CDFs for each random
            %    variable's realization and numpoints is the number
            %    of random vectors. Entries must be in [0,1].
            % N: The number of grid points to use for dividing each axis
            %    of the unit square. The resulting grid of points will be
            %    of size N-by-N. Defaults to 100.
            %
            %OUTPUT:
            % jdist: A 1-by-N^2 vector of joint distribution values.
            % points: A 2-by-N^2 vector of points corresponding to the
            %         probabilities in jprob.
            %
            %EXAMPLE 1: Generates a linearly correlated vector of Gaussian
            %           marginals and plots the copula distribution.
            % % Generate some data
            % numd1 = 1e5;
            % data = randn(2,numd1);
            % sigma = chol([1,0.5;0.5,1],'lower');
            % data = sigma*data;
            % d1 = data(1,:);
            % d2 = data(2,:);
            % u1 = GaussianD.CDF(d1);
            % u2 = GaussianD.CDF(d2);
            % [jprob,points] = Empirical2C.CDF([u1;u2],50);
            % X = reshape(points(1,:),50,50);
            % Y = reshape(points(2,:),50,50);
            % surf(X,Y,reshape(jprob,50,50))
            %
            %December 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if(nargin<2||isempty(N))
                % Generate a uniform grid of (0,1)x(0,1), 100x100 points
                N=100;
            end

            x = linspace(0,1,N);
            [X,Y] = meshgrid(x);

            u1 = u(1,:);
            u2 = u(2,:);

            % Compute copula distribution values
            points = [X(:),Y(:)]';
            numpoints = size(points,2);
            counts = zeros(1,numpoints);
            for idx = 1:numpoints
                counts(idx) = sum((u1<points(1,idx))&(u2<points(2,idx)));
            end
            jdist = counts/sum(counts);
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