function a=solveVandermondeSystem(param1,f,isTransposed)
%%SOLVEVANDERMONDESYSTEM Solve The system of real equations V'*a=f or V*a=f
%               where V is a Vandermonde matrix and f is a vector. The
%               reason this function exists, versus just using V'\f or V\f
%               in Matlab, is because Vandermonde matrices are often poorly
%               conditioned. Thus, this function is usually more accurate
%               than one for solving a generic system of linear equations.
%               Systems of the form V'*a=f arise when performing polynomial
%               interpolation. Those of the form V*a=f arise when
%               determining the weights in a quadrature formula when
%               moments are given. The structure of a Vandermonde matrix is
%               described below.
%
%INPUTS: param1 This can either be an (n+1)X(n+1) real Vandermonde matrix
%               V having the format described below, or this can be an
%               (n+1)X1 or 1X(n+1) vector holding the elements of the
%               second row of the Vandermonde matrix. All entries must be
%               real.
%             f The (n=1)X1 vector in the equation V'*a=f or V*a=f. All
%               entries must be real.
%   isTransposed A boolean parameter. If true, then one desires the solution
%               to the system V'*a=f. If false, then one desires the
%               solution to the system V*a=f. If omitted or an empty matrix
%               is passed, the default is true.
%
%OUTPUTS: a The solution to the system V*a=f or V'*A=f, depending on
%           isTransposed.
%
%A Vandermonde matrix V has the form
%V=[1,      1,      1,    ..., 1;
%   x0,     x1,     x2,   ..., xn;
%   x0^2,   x1^2,   x2^2, ..., xn^2;
%   ...     ...     ...   ...  ...
%   x0^n,   x1^n,   x2^n, ..., xn^n];
%
%This function uses Algorithm 4.6.1 of [1] for the system V'*a=f and
%algorithm 4.6.2 of [1] for the system V*a=f.
%
%As an example, consider a Vandermonde matrix with second row:
% x=[0,1e-3,-1e-3,2e-2,-2e-2];
% V=zeros(5,5);
% V(:,1)=x(1).^(0:4);
% V(:,2)=x(2).^(0:4);
% V(:,3)=x(3).^(0:4);
% V(:,4)=x(4).^(0:4);
% V(:,5)=x(5).^(0:4);
% %We want the solution to V*x=f where
% f=[1;2;3;4;5];
% %The error using this algorithm is 
% V*solveVandermondeSystem(x,f,false)-f
% %The error doing it the generic way is
% V*(V\f)-f
%One can see that the error using this algorithm in on the order of 2e-8,
%but is on the order of 1e-5 using matrix inversion. Larger matrices can
%exibit significantly larger differences.
%
%As noted in [1] in Problem 4.6.1 in Chapter 4.6, the determinant of a
%Vandermonde matrix is prod_{n>=i>j>=0}(x(i+1)-x(j+1)). For the above
%example, the determinant is on the order of 5e-21. 
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(param1)-1;
if(size(param1,1)==size(param1,2))
    %Assume a Vandermonde matrix was passed.
    x=param1(2,:);
else
    %Assume the second row of a Vandermonde matrix was passed.
    x=param1;
end

if(nargin<3||isempty(isTransposed))
    isTransposed=true;
end

if(isTransposed==true)
    for k=0:(n-1)
        for i=n:-1:(k+1)
            f(i+1)=(f(i+1)-f(i+1-1))/(x(i+1)-x(i+1-k-1));
        end
    end

    for k=(n-1):-1:0
        for i=k:(n-1)
            f(i+1)=f(i+1)-f(i+1+1)*x(k+1);
        end
    end
else
    for k=0:(n-1)
        for i=n:-1:(k+1)
            f(i+1)=f(i+1)-x(k+1)*f(i+1-1);
        end
    end

    for k=(n-1):-1:0
        for i=(k+1):n
            f(i+1)=f(i+1)/(x(i+1)-x(i+1-k-1));
        end
        
        for i=k:(n-1)
            f(i+1)=f(i+1)-f(i+1+1);
        end
    end
end

a=f;
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
