classdef Gumbel2C
%%GUMBEL2C A collection of functions for the bivariate Gumbel copulae.
%         Implemented functions include: PDF, CDF, tau.
%
%REFERENCES:
%[1] H. Joe and D. Kurowicka, Dependence Modeling: Vine Copula Handbook.
%    World Scientific, 2011.
%[2] C. Brechmann and U. Schepsmeier, "Dependence modeling with c- and
%    d-vine copulas: The r-package cdvine," Jan. 2011.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    methods(Static)
        function jprob = PDF(u,theta,rot)
            %PDF Generates a joint probability for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % theta: A scalar coupling parameter in the interval
            %        [1,inf).
            % rot: A scalar specifying one of four rotations of the
            %      standard bivariate Gumbel copula. Allowed values are 0,
            %      90, 180, and 270. Default is 0.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %EXAMPLE 1: Make a contour plot of the density.
            %x = linspace(0,1,1e2);
            %[X,Y] = meshgrid(x);
            %Z = zeros(length(X),length(Y));
            %for idx = 1:length(X)
            %    for iidx = 1:length(Y)
            %        Z(idx,iidx) = Gumbel2C.PDF([X(idx,iidx);Y(idx,iidx)],3);
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
            % Z(idx,iidx) = Gumbel2C.PDF([X(idx,iidx);Y(idx,iidx)],3);
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
            elseif theta < 1
                error('Theta parameter cannot be less than 1.')
            end
            if ~exist('rot','var') || isempty(rot)
                rot = 0;
            end
            
            if rot==0
                %No modifications necessary.
            elseif rot==90
                u(1,:) = 1-u(1,:);
            elseif rot==180
                u(1,:) = 1-u(1,:);
                u(2,:) = 1-u(2,:);
            elseif rot==270
                u(2,:) = 1-u(2,:);
            else
                error('Invalid rotation. Allowed values are 0, 90, 180, and 270.')
            end
            
            a = -log(u(1,:));
            b = -log(u(2,:));
            
            jprob = Gumbel2C.CDF(u,theta).*(a.*b).^(theta-1)...
                    .*((a.^theta+b.^theta).^(1/theta)+theta-1)...
                    ./(u(1,:).*u(2,:).*(a.^theta+b.^theta).^(2-(1/theta)));
                
            %NaN results from 0*inf in the calculation.
            %These should be 0s.
            nanIdxs = isnan(jprob);
            jprob(nanIdxs) = 0;
            
        end
        
        function jdist = CDF(u,theta,rot)
            %CDF Generates a joint distribution for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % theta: A scalar coupling parameter in the interval
            %        [1,inf).
            % rot: A scalar specifying one of four rotations of the
            %      standard bivariate Gumbel copula. Allowed values are 0,
            %      90, 180, and 270. Default is 0.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint distribution values.
            %
            %EXAMPLE 1: Make a contour plot of the distribution.
            %x = linspace(0,1,1e2);
            %[X,Y] = meshgrid(x);
            %Z = zeros(length(X),length(Y));
            %for idx = 1:length(X)
            %    for iidx = 1:length(Y)
            %        Z(idx,iidx) = Gumbel2C.CDF([X(idx,iidx);Y(idx,iidx)],3);
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
            % Z(idx,iidx) = Gumbel2C.CDF([X(idx,iidx);Y(idx,iidx)],3);
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
            elseif theta < 1
                error('Theta parameter cannot be less than 1.')
            end
            if ~exist('rot','var') || isempty(rot)
                rot = 0;
            end
            
            if rot==0
                %No modifications necessary.
                a = -log(u(1,:));
                b = -log(u(2,:));
            
                jdist = exp(-(a.^theta+b.^theta).^(1/theta));
            elseif rot==90
                a = -log(1-u(1,:));
                b = -log(u(2,:));
            
                jdist = u(2,:)-exp(-(a.^theta+b.^theta).^(1/theta));
            elseif rot==180
                a = -log(1-u(1,:));
                b = -log(1-u(2,:));
            
                jdist = u(1,:)+u(2,:)-1+exp(-(a.^theta+b.^theta).^(1/theta));
            elseif rot==270
                a = -log(u(1,:));
                b = -log(1-u(2,:));
            
                jdist = u(1,:)-exp(-(a.^theta+b.^theta).^(1/theta));
            else
                error('Invalid rotation. Allowed values are 0, 90, 180, and 270.')
            end
            
        end
        
        function t = tau(theta)
            %%TAU Compute Kendall's tau for the copula.
            %
            %INPUTS:
            % theta: A scalar coupling parameter in the interval
            %        [1,inf).
            %
            %OUTPUTS:
            % t: The scalar value of Kendall's tau.
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('theta','var') || isempty(theta)
                theta = 1;
            elseif theta < 1
                error('Theta parameter cannot be less than 1.')
            end
            
            t = 1-1/theta;
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
