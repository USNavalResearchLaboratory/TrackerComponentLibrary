function J2=numDiff2(x,f,fDim,N,epsilon)
%%NUMDIFF2 Numerically compute the second derivatives using central
%          difference formulae based on Lagrange interpolating polynomials.
%          The algorithm is made to alert (and not crash) in the event that
%          f fails. Specifically, if f returns an empty matrix or it any of
%          the components of f are NaNs, then the numDiff function will
%          terminate early, returning an empty matrix to indicate failure.
%          The function only computes second derivatives, not considering
%          cross terms (it does not find a Hessian matrix). Use numDiffHess
%          to get a Hessian matrix with cross terms.
%
%INPUTS: x The xDimX1 vector or scalar point at which the second derivative
%          of the (possibly vector) function is desired.
%        f The scalar or vector function that is to be differentiated. The
%          function f must take x as its parameter and its output should be
%          a scalar or a column vector.
%     fDim The dimensionality of the output of f.
%        N A number >=1 specifying the order of the second derivative
%          approximation. Values for n=1 through 3 are explicitly coded in.
%          For values 3 and above, the coefficients of the second
%          derivative of the Lagrange interpolating polynomial are
%          explicitly solved. If omitted, a value of N=1 is used.
%  epsilon A scalar or xDimX1 vector quantity specifying the finite step
%          size used for numerical differentiation. If a scalar value is
%          given, that value is used for differentiating with respect to
%          elements of xDim. If an xDimX1 value is given, then the
%          corresponding element of epsilon is used to differentiate each
%          element of x. If epsilon is omitted, then
%          epsilon=max(1e-5*x,1e-7); is used.
%
%OUTPUTS: J2 A fDimXxDim matrix of second derivatives. Each column is the
%            derivative vector of f with respect to the corresponding
%            element of x. If at any point the function f returned a NaN or
%            an empty matrix, J2 will be an empty matrix.
%
%The function is similar to the numDiff function. Central-difference
%numerical differentiation is discussed in terms of Lagrange interpolating
%polynomials in Chapter 4.1 of [1]. The Lagrange interpolating polynomials
%themselves are discussed in Chapter 3 of [1]. In general, the coefficients
%of the interpolating polynomial for first derivative central difference
%numerical differentiation come from Equation 4.2, which  expresses them in
%terms of the first derivative of a Lagrange interpolating polynomial. This
%function just extends the concept by using the second-derivative
%coefficients.
%
%REFERENCES:
%[1] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/ Cole, 2011.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
   N=1; 
end

switch(N)
    case 1
        a=[1; -2; 1];
    case 2
        a=[-1/12; 4/3; -5/2; 4/3; -1/12];
    case 3
        a=[1/90; -3/20; 3/2; -49/18; 3/2; -3/20; 1/90];
    otherwise
        [~,~,d2Li]=LagrangeInterpPoly(-N:1:N);
        a=d2Li(end,:);
end

xDim=size(x,1);

if(nargin<5)
    %If epsilon is not specified, then use some ad-hoc default value
    epsilon=max(1e-5*x,1e-7);
end

if(isscalar(epsilon))
   epsilon=repmat(epsilon,[xDim,1]); 
end

J2=zeros(fDim,xDim);
for curEl=1:xDim
    epsCur=epsilon(curEl);
    
    curIdx=1;
    for curP=-N:1:N
        xP=x;
        xP(curEl)=xP(curEl)+curP*epsCur;
        fxP=f(xP);
        if(isempty(fxP)||any(isnan(fxP)))
            J2=[];
            return;
        end
        J2(:,curEl)=J2(:,curEl)+a(curIdx)*fxP;
        curIdx=curIdx+1;
    end
    J2(:,curEl)=J2(:,curEl)/(epsCur^2);
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
