classdef Gaussian2C
%%GAUSSIAN2C A collection of functions for the bivariate Gaussian copulae.
%            The bivariate log-normal and Gaussian copulae have been proven
%            equivalent due to the monotonic relationship between the two
%            distributions [2].
%            Implemented functions include: PDF, CDF, tau, tau2CorrMat.
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
        function jprob = PDF(u,R)
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
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint probabilities.
            %
            %EXAMPLE 1: Make a contour plot of the density.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            %     for iidx = 1:length(Y)
            %         Z(idx,iidx) = Gaussian2C.PDF([X(idx,iidx);Y(idx,iidx)],0.5);
            %     end
            % end
            % contourf(X,Y,Z)
            %
            %EXAMPLE 2: Generate a joint density plot from random samples
            %           with histograms.
            % rvecs = randn(2,1e4);
            % rho = 0.65;
            % R = chol([1,rho;rho,1],'lower');
            % for idx = 1:size(rvecs,2)
            %     rvecs(:,idx) = R*rvecs(:,idx);
            % end
            % c = corrCoeffSpearman(rvecs(1,:),rvecs(2,:)); %diagonal entries close to rho
            % u1 = GaussianD.CDF(rvecs(1,:));
            % u2 = GaussianD.CDF(rvecs(2,:));
            % pdf = Gaussian2C.PDF([u1;u2],rho).*GaussianD.PDF(rvecs(1,:)).*GaussianD.PDF(rvecs(2,:));
            % fig = jointPlot2D(rvecs(1,:),rvecs(2,:),pdf);
            %
            %EXAMPLE 3: Generate an exact joint density plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            % for iidx = 1:length(Y)
            % Z(idx,iidx) = Gaussian2C.PDF([X(idx,iidx);Y(idx,iidx)],0.95);
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
            
            a = sqrt(2)*erfinv(2*u(1,:)-1);
            b = sqrt(2)*erfinv(2*u(2,:)-1);
            c = (0.5*(2*a.*b*rho-rho^2*(a.^2+b.^2))/(1-rho^2));
            
            jprob = (1/sqrt(1-rho^2))*exp(c);
            
            %NaN results from subtraction of infinities in the exponential.
            %These should be 0s.
            nanIdxs = isnan(jprob);
            jprob(nanIdxs) = 0;
            
        end
        
        function jdist = CDF(u,R)
            %CDF Generates a joint distribution for column vector u or
            %    a row vector of joint probabilities for matrix u
            %    with each column as the input vector.
            %
            %INPUTS: 
            % u: A 2-by-n matrix where each column represents a
            %    2-dimensional random vector and n is the number of
            %    random vectors. Entries must be in [0,1].
            % R: A 2-by-2 correlation matrix, or a scalar which is the
            %    correlation.
            %
            %OUTPUTS: jprob: A 1-by-n vector of joint distribution values.
            %
            %Note: This function uses a Monte Carlo method for
            %      approximating the CDF. It will take a few seconds to run
            %      the examples below.
            %
            %EXAMPLE 1: Make a contour plot of the distribution.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            %     for iidx = 1:length(Y)
            %         Z(idx,iidx) = Gaussian2C.CDF([X(idx,iidx);Y(idx,iidx)],0.65);
            %     end
            % end
            % contourf(X,Y,Z)
            %
            %EXAMPLE 2: Generate an exact joint distribution plot.
            % x = linspace(0,1,1e2);
            % [X,Y] = meshgrid(x);
            % Z = zeros(length(X),length(Y));
            % for idx = 1:length(X)
            %     for iidx = 1:length(Y)
            %         Z(idx,iidx) = Gaussian2C.CDF([X(idx,iidx);Y(idx,iidx)],0.65);
            %         X(idx,iidx) = GumbelD.invCDF(X(idx,iidx),1,2);
            %         Y(idx,iidx) = GaussianD.invCDF(Y(idx,iidx));
            %         Z(idx,iidx) = Z(idx,iidx);
            %     end
            % end
            % tiledlayout(1,2)
            % nexttile
            % surf(X,Y,Z)
            % Zsmooth = zeros(size(Z));
            % Zsmooth(:,end) = Z(:,end);
            % Zsmooth(:,1) = Z(:,1);
            % Zsmooth(end,:) = Z(end,:);
            % Zsmooth(1,:) = Z(1,:);
            % for idx = 2:size(Z,1)-1
            %     for iidx = 2:size(Z,2)-1
            %         Zsmooth(idx,iidx) = sum(Z(idx-1:idx+1,iidx-1:iidx+1),'all')/9;
            %     end
            % end
            % nexttile
            % surf(X,Y,Zsmooth)
            %
            %October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if ~exist('u','var') || isempty(u)
                jdist = [];
                return
            end
            if ~exist('R','var') || isempty(R)
                rho = 0;
            elseif ~isscalar(R)
                rho = R(2,1);
            else
                rho = R;
            end
             
            R = [1,rho;rho,1]^2;
            
            a = sqrt(2)*erfinv(2*u(1,:)-1);
            b = sqrt(2)*erfinv(2*u(2,:)-1);

            %Though this approximation could be used, the function
            %bivarGaussRectangleCDF produces a more exact result. Moreover,
            %it can handle -Inf lower bound, whereas the approximation
            %would have to set some arbitrary llim as the lower bound.
            %  jdist = arrayfun(@(ele1,ele2)GaussianD.integralOverRegion([0;0],R,[llim;llim],[ele1;ele2]),a,b);
            jdist = arrayfun(@(ele1,ele2)bivarGaussRectangleCDF([-Inf;-Inf],[ele1;ele2],[0;0],R),a,b);
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
            %EXAMPLE: Shows Gaussian distribution of independent
            %         correlation parameter estimates using Kendall's
            %         tau.
            % for idx = 1:1e3
            %    x = GaussianD.rand(1e3,[0;0],[1,0.3;0.3,1]);
            %    u = GaussianD.CDF(x);
            %    s1 = GammaD.invCDF(u(1,:),2,1);
            %    s2 = GammaD.invCDF(u(2,:),4,4);
            %    tauEst(idx) = corr(s1',s2','Type','Kendall');
            % end
            % histogram(tauEst)
            % mu = mean(tauEst);
            % sigma = std(tauEst);
            % xline(mu,'k','LineWidth',2,'Label','\mu');
            % xline(mu-sigma,'r--','LineWidth',2,'Label','\mu-\sigma');
            % xline(mu+sigma,'r--','LineWidth',2,'Label','\mu+\sigma');
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
        
        function R = tau2corrMat(t)
            %%TAU2CORRMAT Compute a 2-by-2 correlation matrix estimate
            %             from Kendall's tau.
            %
            %INPUTS:
            % t: Kendall's tau
            %
            %OUTPUTS:
            % R: A 2-by-2 correlation matrix
            %
            %EXAMPLE 1: Generates Gamma distributed marginals correlated
            %           via a Gaussian copula and estimates the correlation
            %           matrix.
            % x = GaussianD.rand(1e4,[0;0],[1,0.3;0.3,1]);
            % u = GaussianD.CDF(x);
            % s1 = GammaD.invCDF(u(1,:),2,1);
            % s2 = GammaD.invCDF(u(2,:),4,4);
            % tau = corr(s1',s2','Type','Kendall');
            % Gaussian2C.tau2corrMat(tau)
            %
            %EXAMPLE 2: Shows Gaussian distribution of independent
            %           correlation parameter estimates using Kendall's
            %           tau.
            % for idx = 1:1e3
            %    x = GaussianD.rand(1e3,[0;0],[1,0.3;0.3,1]);
            %    u = GaussianD.CDF(x);
            %    s1 = GammaD.invCDF(u(1,:),2,1);
            %    s2 = GammaD.invCDF(u(2,:),4,4);
            %    tauEst(idx) = corr(s1',s2','Type','Kendall');
            %    rho = Gaussian2C.tau2corrMat(tauEst(idx));
            %    theta(idx) = rho(1,2);
            % end
            % histogram(theta)
            % mu = mean(theta);
            % sigma = std(theta);
            % xline(mu,'k','LineWidth',2,'Label','\mu');
            % xline(mu-sigma,'r--','LineWidth',2,'Label','\mu-\sigma');
            % xline(mu+sigma,'r--','LineWidth',2,'Label','\mu+\sigma');
            %
            %April 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            rho = sin(t*pi/2);
            R = [1,rho;rho,1];
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
