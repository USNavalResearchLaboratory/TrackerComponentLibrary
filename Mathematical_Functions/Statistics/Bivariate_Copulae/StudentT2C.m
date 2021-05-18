classdef StudentT2C
%%STUDENTT2C A collection of functions for the bivariate Student-t copulae.
%           Implemented functions include: PDF, tau.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function jprob = PDF(u,R,nu)
            %PDF Generates a joint probability for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % R: A 2-by-2 correlation matrix, or a scalar which is the
            %    correlation.
            % nu: Scalar degrees of freedom.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %EXAMPLE 1: Make a contour plot of the density.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            %     for iidx = 1:length(Y)
            %         Z(idx,iidx) = StudentT2C.PDF([X(idx,iidx);Y(idx,iidx)],0.5,3);
            %     end
            % end
            % contourf(X,Y,Z)
            %
            %EXAMPLE 2: Generate an exact joint density plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            % for iidx = 1:length(Y)
            % Z(idx,iidx) = StudentT2C.PDF([X(idx,iidx);Y(idx,iidx)],0.95,3);
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
            if ~exist('R','var') || isempty(R)
                rho = 0;
            elseif ~isscalar(R)
                rho = R(2,1);
            else
                rho = R;
            end
            
            a = StudentTD.invCDF(u(1,:),nu);
            b = StudentTD.invCDF(u(2,:),nu);
            c = (a.^2+b.^2-2*rho*a.*b)/(nu*(1-rho^2));
            
            jprob = (gamma(0.5*nu+1)/(gamma(nu/2)*nu*pi*sqrt(1-rho^2)))...
                    ./(StudentTD.PDF(a,0,1,nu).*StudentTD.PDF(b,0,1,nu).*(1+c).^(0.5*nu+1));
                
            %NaN results from 0*inf in the calculation.
            %These should be 0s.
            nanIdxs = isnan(jprob);
            jprob(nanIdxs) = 0;
            
        end
        
        function t = tau(R)
            %%TAU Compute Kendall's tau for the copula.
            %
            %INPUTS:
            % R: A 2-by-2 correlation matrix, or a scalar which is the
            %    correlation.
            %
            %OUTPUTS:
            % t: The scalar value of Kendall's tau.
            %
            %Note: This does not depend on the degrees of freedom nu.
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('R','var') || isempty(R)
                rho = 0;
            elseif ~isscalar(R)
                rho = R(2,1);
            else
                rho = R;
            end
            
            t = 2*arcsin(rho)/pi;
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
