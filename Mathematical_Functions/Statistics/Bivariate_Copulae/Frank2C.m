classdef Frank2C
%%FRANK2C A collection of functions for the bivariate Frank copulae.
%         Implemented functions include: PDF, CDF, tau.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function jprob = PDF(u,theta)
            %PDF Generates a joint probability for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % theta: A scalar coupling parameter in the interval
            %        (-inf,inf)\{0}.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %EXAMPLE 1: Make a contour plot of the density.
            %x = linspace(0,1,1e2);
            %[X,Y] = meshgrid(x);
            %Z = zeros(length(X),length(Y));
            %for idx = 1:length(X)
            %    for iidx = 1:length(Y)
            %        Z(idx,iidx) = Frank2C.PDF([X(idx,iidx);Y(idx,iidx)],0.5);
            %    end
            %end
            %surf(X,Y,Z)
            %
            %EXAMPLE 2: Generate an exact joint density plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            % for iidx = 1:length(Y)
            % Z(idx,iidx) = Frank2C.PDF([X(idx,iidx);Y(idx,iidx)],0.95);
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
            if ~exist('theta','var') || isempty(theta)
                theta = 1;
            elseif theta == 0
                error('Theta parameter cannot equal 0.')
            end
            
            a = exp(-theta*u(1,:));
            b = exp(-theta*u(2,:));
            c = exp(-theta);
            
            jprob = theta*a.*b.*(1-c)./(c-a-b+a.*b).^2;
            
        end
        
        function jdist = CDF(u,theta)
            %CDF Generates a joint distribution value for column vector u or
            %    a row vector of joint distribution values for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % theta: A scalar coupling parameter in the interval
            %        (-inf,inf)\{0}.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %EXAMPLE 1: Make a contour plot of the density.
            %x = linspace(0,1,1e2);
            %[X,Y] = meshgrid(x);
            %Z = zeros(length(X),length(Y));
            %for idx = 1:length(X)
            %    for iidx = 1:length(Y)
            %        Z(idx,iidx) = Frank2C.PDF([X(idx,iidx);Y(idx,iidx)],0.5);
            %    end
            %end
            %surf(X,Y,Z)
            %
            %EXAMPLE 2: Generate an exact joint distribution plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            % for iidx = 1:length(Y)
            % Z(idx,iidx) = Frank2C.CDF([X(idx,iidx);Y(idx,iidx)],0.95);
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
            if ~exist('theta','var') || isempty(theta)
                theta = 1;
            elseif theta == 0
                error('Theta parameter cannot equal 0.')
            end
            
            a = exp(-theta*u(1,:));
            b = exp(-theta*u(2,:));
            c = exp(-theta);
            
            jdist = -(1/theta)*log(1+(a-1).*(b-1)/(c-1));
            
        end
        
        function t = tau(theta)
            %%TAU Compute Kendall's tau for the copula.
            %
            %INPUTS:
            % theta: A scalar coupling parameter in the interval
            %        (-inf,inf)\{0}.
            %
            %OUTPUTS:
            % t: The scalar value of Kendall's tau.
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('theta','var') || isempty(theta)
                theta = 1;
            elseif theta == 0
                error('Theta parameter cannot equal 0.')
            end
            
            t = 1+4*(Debye(theta)-1)/theta;
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
