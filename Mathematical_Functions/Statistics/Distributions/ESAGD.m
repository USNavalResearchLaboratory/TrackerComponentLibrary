classdef ESAGD
    %%ESAGD The ellipically symmetric angular Gaussian distribution (ESAG)
    %       is asymmetric projection of a multivariate Gaussian onto the
    %       unit sphere. This class provides functions for computations
    %       relevant to the ESAG in 3 dimensions. The class has two
    %       parameter vectors, mu and gamma. Mu is the mean of the
    %       multivariate Gaussian being projected. Gamma is a vector of
    %       two quantities, first the ratio of nonunit eigenvalues and
    %       second the rotation of the principle axes defined by the
    %       tangent plane at the projection of mu onto the unit sphere.
    %
    %Implemented methods are: getGamma, PDF, derivedParams, rand, fit
    %
    %REFERENCES:
    %[1] P. Paine, S. P. Preston, M. Tsagris, and A. T. Wood, "An
    %    elliptically symmetric angular gaussian distribution," Statistics
    %    and Computing, vol. 28, no. 3, pp. 689-697, 2018.
    %
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    methods(Static)
        function gamma = getGamma(mu,rho,psi)
            %%GETGAMMA Computes the gamma parameter as defined in Eq. (17)
            %          of [1], though there is a difference between the psi
            %          as defined in [1] and as used here. See description
            %          below.
            %
            %INPUTS:
            % mu: A 3-by-1 vector specifying the mean of the density.
            % rho: The smallest eigenvalue of V. This should be in the
            %      interval (0,1].
            % psi: A scalar rotation in the interval [0,pi] which describes
            %      how the az/el axes are rotated about the radial axis.
            %      This should be in radians.
            %
            %      Note this is not the same as psi in [1]. In [1], a psi
            %      of 0 does not mean 0 rotation about the radial axis if
            %      mu is off the plane defined by the second and third
            %      Cartesian coordinates (mu(1) nonzero). The conversion
            %      needed to offset this induced rotation of the axes
            %      tangential to the sphere's surface is an addition of the
            %      angle rotated in the plane defined by the first and
            %      third Cartesian axes.
            %
            %OUTPUTS:
            % gamma: A 2-by-1 vector containing the gamma parameters as
            %        defined in Eq. (17) of [1].
            %
            %EXAMPLE: Plots an ESAG distribution on the unit sphere using a
            %         given mu and gamma computed with this function and
            %         some choices for rho and psi.
            % mu = [3;3;5];
            % gam = ESAGD.getGamma(mu,0.5,pi/2);
            % [X,Y,Z] = sphere(100);
            % Cmeas = reshape(ESAGD.PDF([X(:),Y(:),Z(:)]',mu,gam),size(X,1),size(Y,1));
            % surf(X,Y,Z,'CData',Cmeas,"EdgeColor","none")
            % view([110,50])
            % axis('equal')
            % title("Angular Measurement Projected onto Unit Sphere")
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            aor = atan(mu(1)/mu(2));
            if mu(3)>0
                psi = psi+aor;
            else
                psi = psi-aor;
            end
            gamma = [(1/rho-rho)*cos(2*psi)/2;
                (1/rho-rho)*sin(2*psi)/2];
        end
        
        function pdf = PDF(x,mu,gamma)
            %PDF Computes the probability density of an ESAG distribution
            %    with parameters mu and gamma.
            %
            %INPUTS:
            % x: A xDim-by-N matrix containing N vectors at which the
            %    density will be evaluated.
            % mu: A 3-by-1 vector specifying the mean of the density.
            % gamma: A 2-by-1 vector specifying the scaling and rotation of
            %        the density (see getGamma method if the scaling and
            %        rotation angles are known). gamma(1) should be in
            %        (0,1] and gamma(2) should be in radians and in [0,pi].
            %
            %OUTPUTS:
            % pdf: A N-by-1 vector of density values corresponding to the
            %      column vectors in x.
            %
            %EXAMPLE: Plots an ESAG distribution on the unit sphere.
            % mu = [3;3;5];
            % gam = ESAGD.getGamma(mu,0.5,pi/2);
            % [X,Y,Z] = sphere(100);
            % Cmeas = reshape(ESAGD.PDF([X(:),Y(:),Z(:)]',mu,gam),size(X,1),size(Y,1));
            % surf(X,Y,Z,'CData',Cmeas,"EdgeColor","none")
            % view([110,50])
            % axis('equal')
            % title("Angular Measurement Projected onto Unit Sphere")
            %
            %This function is a modified version of one provided by
            %S. P. Preston in connection with reference [1]. The
            %original function can be found on that author's website.
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            pdf = PDF(x,mu,gamma);
        end
        
        function [invV, ksi1, ksi2, ksi3, mu0] = derivedParams(mu,gamma)
            %%DERIVEDPARAMS Computes the various derived parameters defined
            %               in [1].
            %
            %INPUTS:
            % mu: A 3-by-1 vector specifying the mean of the density.
            % gamma: A 2-by-1 vector specifying the scaling and rotation of
            %        the density (see getGamma method if the scaling and
            %        rotation angles are known). Defaults to [1;0], the
            %        standard isometric angular Gaussian (IAG) about mu.
            %        gamma(1) should be in (0,1] and gamma(2) should be in
            %        radians and in [0,pi].
            %
            %OUTPUTS:
            % invV: The 3-by-3 inverse of the covariance matrix V as
            %       defined in Eq. (18).
            % ksi1: The 3-by-1 principle axis rotated from the elevation
            %       axis.
            % ksi2: The 3-by-1 principle axis rotated from the azimuth
            %       axis.
            % ksi3: The 3-by-1 axis normal to the sphere and parallel with
            %       the vector mu.
            % mu0: The magnitude of the mu vector restricted to the second
            %      and third coordinates.
            %
            %This function is a modified version of one provided by
            %S. P. Preston in connection with reference [1]. The
            %original function can be found on that author's website.
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            
            % calls the private function defined below
            [invV, ksi1, ksi2, ksi3, mu0] = derivedParams(mu,gamma);
        end
        
        function points = rand(N,mu,gamma)
            %%RAND Generates a collection of random samples from the ESAG
            %      described by parameters mu and gamma.
            %
            %INPUTS:
            % N: A scalar number of points to be sampled.
            % mu: A 3-by-1 vector specifying the mean of the density.
            % gamma: A 2-by-1 vector specifying the scaling and rotation of
            %        the density (see getGamma method if the scaling and
            %        rotation angles are known).
            %
            %OUTPUTS:
            % points: A 3-by-N matrix of points sampled from the
            %         distribution parameterized by mu and gamma.
            %
            %This function is a modified version of one provided by
            %S. P. Preston in connection with reference [1]. The
            %original function can be found on that author's website.
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            invV = derivedParams(mu,gamma);
            x = GaussianD.rand(N,mu,invV\eye(3));
            points = x./vecnorm(x);
        end
        
        function [Bhat, optimOpts] = fit(Y,X,B0,params2Fit,optimOpts)
            %%FIT Fits ESAG parameters [mu;gamma] to the given data (X,Y).
            %
            %INPUTS:
            % Y: A 3-by-N matrix of sample vectors.
            % X: A p-by-N matrix of input parameters for each corresponding
            %    sample vector in Y. Defaults to a matrix of ones.
            % B0: A 5-by-1 vector of initialization parameter estimates.
            %     Defaults to a vector of ones.
            % params2Fit: A 5-by-1 boolean vector indicating which
            %             parameters to optimize over (with 1) and which to
            %             leave as fixed (with 0). The order is [mu;gamma].
            %             Defaults to a vector of ones.
            % optimOpts: A collection of arguments to pass to the
            %            fminsearch function to modify its behavior.
            %            Defaults are set with the following: 
            %            optimset('Display','off',...
            %                     'MaxFunEvals',100000,...
            %                     'MaxIter',100000,...
            %                     'TolX',1e-8,...
            %                     'TolFun',1e-8)
            %
            %OUTPUTS:
            % Bhat: A 5-by-1 vector of the final estimated parameters.
            % optimOpts: The optimset parameters used for modifying the
            %            behavior of fminsearch.
            %
            %This function is a modified version of one provided by
            %S. P. Preston in connection with reference [1]. The
            %original function can be found on that author's website.
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            if nargin<2 || isempty(X)
                X = ones(5,size(Y,2));
            end
            if nargin<3 || isempty(B0)
                B0 = ones(5,1);
            end
            if nargin<4
                optimOpts = optimset('Display','off','MaxFunEvals',...
                100000,'MaxIter',100000,'TolX',1e-8,'TolFun',1e-8);
            end
            a0 = B0(parsToFit==1);
            b = B0(parsToFit==0);
            [ahat,optimOpts] = fminsearch(@(a) f(a, b, Y, X, parsToFit),a0,optimOpts);
            Bhat = NaN(5,p);
            Bhat(params2Fit==1) = ahat;
            Bhat(params2Fit==0) = b;
        end
    end
end

function [ksi1,ksi2,ksi3,mu0] = getProjAxes(mu)
    %%GETPROJAXES Computes the unrotated versions of ksi1, ksi2 as
    %             given in Eq. (14) of [1].
    %
    %INPUTS:
    % mu: A 3-by-1 vector specifying the mean of the density.
    %
    %OUTPUTS:
    % ksi1: The first 3-by-1 principle axis.
    % ksi2: The second 3-by-1 principle axis.
    % ksi3: The 3-by-1 axis normal to the sphere and parallel with the
    %       vector mu.
    % mu0: The magnitude of the mu vector restricted to the second and
    %      third coordinates.
    %
    %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
    %
    mu0 = sqrt(mu(2)^2+mu(3)^2);
    ksi1 = [-mu0^2;mu(1)*mu(2);mu(1)*mu(3)]/(mu0*norm(mu));
    ksi2 = [0;-mu(3);mu(2)]/mu0;
    ksi3 = mu/norm(mu);
end

%This functionality is given as a method above, but defined here for use
%in the other class methods.
function [invV, ksi1t, ksi2t, ksi3, mu0] = derivedParams(mu,gamma)
%%DERIVEDPARAMS Computes the various derived parameters defined
%               in [1].
%
%INPUTS:
% mu: A 3-by-1 vector specifying the mean of the density.
% gamma: A 2-by-1 vector specifying the scaling and rotation of
%        the density (see getGamma method if the scaling and
%        rotation angles are known). Defaults to [1;0], the
%        standard isometric angular Gaussian (IAG) about mu.
%        gamma(1) should be in (0,1] and gamma(2) should be in
%        radians and in [0,pi].
%
%OUTPUTS:
% invV: The 3-by-3 inverse of the covariance matrix V as
%       defined in Eq. (18).
% ksi1: The 3-by-1 principle axis rotated from the elevation
%       axis.
% ksi2: The 3-by-1 principle axis rotated from the azimuth
%       axis.
% ksi3: The 3-by-1 axis normal to the sphere and parallel with the
%       vector mu.
% mu0: The magnitude of the mu vector restricted to the second and
%      third coordinates.
%
%This function is a modified version of one provided by
%S. P. Preston in connection with reference [1]. The
%original function can be found on that author's website.
%
%November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%
if nargin<2
    gamma = [1;0];
end
[ksi1,ksi2,ksi3,mu0] = getProjAxes(mu);
ksi1t = ksi1*cos(gamma(2))+ksi2*sin(gamma(2));
ksi2t = -ksi1*sin(gamma(2))+ksi2*cos(gamma(2));
invV = eye(3)+gamma(1)*(ksi1*ksi1'-ksi2*ksi2')+...
    gamma(2)*(ksi1*ksi2'+ksi2*ksi1')+...
    (sqrt(1+gamma(1)^2+gamma(2)^2)-1)*(ksi1*ksi1'+ksi2*ksi2');
end

%This functionality is given as a method above, but defined here for use
%in the other class methods.
function pdf = PDF(x,mu,gamma)
%PDF Computes the probability density of an ESAG distribution
%    with parameters mu and gamma.
%
%INPUTS:
% x: A xDim-by-N matrix containing N vectors at which the
%    density will be evaluated.
% mu: A 3-by-1 vector specifying the mean of the density.
% gamma: A 2-by-1 vector specifying the scaling and rotation of
%        the density (see getGamma method if the scaling and
%        rotation angles are known). gamma(1) should be in
%        (0,1] and gamma(2) should be in radians and in [0,pi].
%
%OUTPUTS:
% pdf: A N-by-1 vector of density values corresponding to the
%      column vectors in x.
%
%EXAMPLE: Plots an ESAG distribution on the unit sphere.
% mu = [3;3;5];
% gam = ESAGD.getGamma(mu,0.5,pi/2);
% [X,Y,Z] = sphere(100);
% Cmeas = reshape(ESAGD.PDF([X(:),Y(:),Z(:)]',mu,gam),size(X,1),size(Y,1));
% surf(X,Y,Z,'CData',Cmeas,"EdgeColor","none")
% view([110,50])
% axis('equal')
% title("Angular Measurement Projected onto Unit Sphere")
%
%This function is a modified version of one provided by
%S. P. Preston in connection with reference [1]. The
%original function can be found on that author's website.
%
%November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%
[~,numSamp] = size(x);
invV = derivedParams(mu,gamma);
pdf = zeros(numSamp,1);
c3 = mu'*mu;
for i = 1:numSamp
    y = x(:,i)/norm(x(:,i));
    c1 = y'*invV*y;
    if c1>=0
        c2 = y'*mu;
        alpha = c2/sqrt(c1);
        M2 = (1+alpha.^2).*GaussianD.CDF(alpha)+alpha.*GaussianD.PDF(alpha);
        pdf(i) = (M2/(2*pi))*c1^(-3/2)*exp(1/2*(alpha^2-c3));
    else
        pdf(i) = -inf;
    end
end
end

function out = f(a, b, Y, X, params2Fit)
%%F A function for optimizing the log-likelihood of ESAG in the fit
%   method.
%
%INPUTS:
% a: 5-by-1 vector of fit parameters.
% b: 5-by-1 vector of fixed parameters.
% Y: A 3-by-N matrix of sample vectors.
% X: A p-by-N matrix of input parameters for each corresponding
%    sample vector in Y.
% params2Fit: A 5-by-1 boolean vector indicating true for any
%             parameters which should be fit.
%
%OUTPUTS:
% out: The log-likelihood of an ESAG variable with the chosen
%      parameters.
%
%This function is a modified version of one provided by
%S. P. Preston in connection with reference [1]. The
%original function can be found on that author's website.
%
%November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%
[p,~] = size(X);
[~,n] = size(Y);

B = zeros(5,p);
B(params2Fit==1) = a;
B(params2Fit==0) = b;

LL = 0;

for i = 1:n
   mu = B(1:3,:)*X(:,i);
   gam = B(4:5,:)*X(:,i);
   param = [mu;gam];
   LL = LL + log(PDF(Y(:,i),param));
end
out = -LL;
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