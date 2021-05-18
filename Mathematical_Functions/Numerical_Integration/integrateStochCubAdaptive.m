function [I,V,N]=integrateStochCubAdaptive(func,mu,SR,algorithm,epsVal2,NMax)
%%INTEGRATESTOCHCUBADAPTIVE Perform approximate Monte Carlo numerical
%                integration of a function times a Gaussian PDF using a
%                method that can be compared to importance sampling in that
%                random sets of approximate cubature points are generated
%                rather than more traditional Monte Carlo appraoches.
%
%INPUTS: func The handle to the function that is multiplies by a Gaussian
%             PDF before integration is performed. The function can take a
%             multidimensional input, but must produce a real univariate
%             output.
%          mu The numDimX1 mean of the Gaussian PDF.
%          SR The numDimXnumDim lower-triangular square root of the
%             covariance matrix of the Gaussian PDF.
%   algorithm An optional parameter specifying the algorithm to use.
%             Possible values are:
%             0 Use the first-order method of [1].
%             1 Use the third-order method of [1].
%             2 (The default if omitted or an empty matrix is passed) Use
%               the fifth-order method of [1].
%     epsVal2 An optional convergence bound based on the approximate
%             variance of the estimate. The default if this parameter is
%             omitted or an empty matrix is passed is 1e-2.
%        NMax The maximum number of iterations to perform. The default
%             value if this parameter is omitted or an empty matrix is
%             passed is 10e3.
%
%OUTPUTS: I The approximate scalar value of the integral.
%         V The approximate variance of the estimate.
%         N The number of iterations performed.
%
%All three algorithms arise from similar derivations in [1].
%
%EXAMPLE:
%Here, we choose a function whose true integral with a normal weighting
%function can be easily found. In this instance, we just choose a bivariate
%polynomial. 
% f=@(x)(x(2)^4*(x(1)+1)^3-(x(1)+2)^2-(x(2)-4)^3);
% mu=[0;1/2];
% R=[2,-1;
%   -1, 2];
% SR=chol(R,'lower');
% algorithm=2;
% %The exact solution to the problem is known to be 
% exactSol=1917/16
% [I,V,N]=integrateStochCubAdaptive(f,mu,SR,algorithm)
%One will find the result is relatively close, but not equal to the exact
%solution.
%
%REFERENCES:
%[1] A. Genz and J. Monahan, "Stochastic integration rules for infinite
%    regions," SIAM Journal on Scientific Computing, vol. 19, no. 2, pp.
%    426-439, Mar. 1998.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(NMax))
   NMax=10e3; 
end

if(nargin<5||isempty(epsVal2))
   epsVal2=1e-2; 
end

if(nargin<4||isempty(algorithm))
    algorithm=2;
end

numDim=length(mu);

%Deal with the distribution not being a standard normal distribution.
f=@(x)func(SR*x+mu);

switch(algorithm)
    case 0 %Degree 1 spherical-radial integration rule 
        %Step 2
        N=0;
        I=0;
        V=0;
        while(1)
           %Step 3a
           N=N+1;

           %Step3b
           x=randn(numDim,1);

           %Step 3c
           SR=(f(-x)+f(x))/2;
           D=(SR-I)/N;
           I=I+D;
           V=(N-2)*V/N+D^2;

           if(V<epsVal2||N==NMax)
               break;
           end
        end

        %V is approximate sample variance....
    case 1 %Degree 3 spherical-radial integration rule
        %Step 2
        N=0;
        I=0;
        V=0;
        F0=f(zeros(numDim,1));
        e=zeros(numDim,1);
        while(1)
            %Step 3a
            N=N+1;
            SR=0;
            %Step 3b
            Q=randOrthoMat(numDim);
            %Step 3c
            
            rho=ChiD.rand(1,numDim+2);

            %Step 3d
            for j=1:numDim
                e(j)=1;
                x=rho*Q*e;

                SR=SR+f(-x)+f(x);
                e(j)=0;
            end

            %Step 3e
            SR=F0*(1-numDim/rho^2)+SR/(2*rho^2);
            D=(SR-I)/N;
            I=I+D;
            V=(N-2)*V/N+D^2;

            if(V<epsVal2||N==NMax)
                break;
            end
        end
    case 2%Degree 5 spherical-radial integration rule
        m=numDim;
        %Step 2
        N=0;
        I=0;
        V=0;
        F0=f(zeros(m,1));
                
        v=regularNSimplexCoords(m);
        while(1)
            %Step 3a
            N=N+1;
            %Step 3b
            Q=randOrthoMat(m);
            vTilde=Q*v;
            %Step 3c
            r=ChiD.rand(1,2*m+7);
            q=BetaD.rand(1,m+2,3/2);
            rho=r*sin(asin(q)/2);
            delta=r*cos(asin(q)/2);
            Fv=0;
            Fy=0;

            %Step 3d
            for j=1:(m+1)
                x1=rho*vTilde(:,j);
                x2=delta*vTilde(:,j);

                Fv=Fv+(m+2-delta^2)*(f(-x1)+f(x1))/(rho^2*(rho^2-delta^2))...
                     +(m+2-rho^2)*(f(-x2)+f(x2))/(delta^2*(delta^2-rho^2));

                for i=1:(j-1)
                    y=(vTilde(:,j)+vTilde(:,i))/norm(vTilde(:,j)+vTilde(:,i));

                    x1=rho*y;
                    x2=delta*y;

                    Fy=Fy+(m+2-delta^2)*(f(-x1)+f(x1))/(rho^2*(rho^2-delta^2))...
                         +(m+2-rho^2)*(f(-x2)+f(x2))/(delta^2*(delta^2-rho^2));
                end
            end

            %Step 3e
            SR=F0*(1-m*(rho^2+delta^2-(m+2))/(rho^2*delta^2))+(Fv*(7-m)*m^2+4*Fy*(m-1)^2)/(2*(m+1)^2*(m+2));
            D=(SR-I)/N;
            I=I+D;
            V=(N-2)*V/N+D^2;

            if(V<epsVal2||N==NMax)
                break;
            end
        end
    otherwise
        error('Unknown algorithm specified.')
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
