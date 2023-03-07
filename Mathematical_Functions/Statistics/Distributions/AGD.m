classdef AGD
    %%AGD A class containing functions related to the angular Gaussian
    %     distribution (AGD). This is the distribution characterized by
    %     the angular projection of a Gaussian ellipsoid onto the unit
    %     sphere.
    %
    %Implemented methods are: PDF, rand
    %
    %REFERENCES:
    %[1] P. Paine, S. P. Preston, M. Tsagris, and A. T. Wood, "An
    %    elliptically symmetric angular gaussian distribution," Statistics
    %    and Computing, vol. 28, no. 3, pp. 689-697, 2018.
    %
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    methods(Static)
        function pdf = PDF(x,mu,V)
            %PDF Computes the probability density of an AG distribution
            %    with mean and covariance information.
            %
            %INPUTS:
            % x: A xDim-by-N matrix containing N vectors at which the
            %    density will be evaluated.
            % mu: A 3-by-1 vector specifying the mean of the density.
            % V: A 3-by-3 matrix specifying the covariance of the
            %    distribution.
            %
            %OUTPUTS:
            % pdf: A N-by-1 vector of density values corresponding to the
            %      column vectors in x.
            %
            %EXAMPLE: Plots a multivariate Gaussian ellipsoid and its
            %         angular projection on the unit sphere.
            % mu = [3;3;5];
            % spheremu = mu/norm(mu);
            % V = [1,0.3,0;0.2,0.5,0.7;0,0,5]; V = V'*V;
            % 
            % [X,Y,Z] = sphere(100);
            % Ctrack = reshape(AGD.PDF([X(:),Y(:),Z(:)]',mu,V),size(X,1),size(Y,1));
            % 
            % surf(X,Y,Z,'CData',Ctrack,"EdgeColor","none")
            % hold on
            % drawEllipse(mu,V\eye(3),ChiSquareD.invCDF(0.2,3),'FaceAlpha',0.2,'EdgeAlpha',0.2)
            % plot3([mu(1);spheremu(1)],[mu(2);spheremu(2)],[mu(3);spheremu(3)],'-r','LineWidth',2)
            % hold off
            % view([110,50])
            % axis('equal')
            % title("Angular Projection of 3D Gaussian onto Unit Sphere")
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            [xDim,numSamp] = size(x);
            
            invV = V\eye(xDim);
            y = x./repmat(vecnorm(x),xDim,1);
            c1 = zeros(numSamp,1);
            c2 = zeros(numSamp,1);
            for i = 1:numSamp
                c1(i) = y(:,i)'*invV*y(:,i);
                c2(i) = y(:,i)'*invV*mu;
            end
            c3 = mu'*invV*mu;
            alpha = c2./sqrt(c1);
            M2 = (1+alpha.^2).*GaussianD.CDF(alpha)+alpha.*GaussianD.PDF(alpha);
            pdf = (M2/(2*pi)).*c1.^(-3/2).*exp(1/2*(alpha.^2 - c3));
        end
        
        function points = rand(N,mu,V)
            %%RANDFROMV Generates a collection of random samples from the
            %           AGD described by the mean and covariance.
            %
            %INPUTS:
            % N: A scalar number of points to be sampled.
            % mu: A 3-by-1 vector specifying the mean of the density.
            % V: A 3-by-3 matrix specifying the covariance of the
            %    distribution.
            %
            %OUTPUTS:
            % points: A 3-by-N matrix of points sampled from the
            %         distribution parameterized by mu and gamma.
            %
            %November 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
            %
            x = GaussianD.rand(N,mu,V);
            points = x./vecnorm(x);
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