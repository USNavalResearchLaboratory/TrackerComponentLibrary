function [intT,T]=ChebyshevPolyInt(tau,n,tauStart,tauEnd)
%%CHEBYSHEVPOLYDERIV Compute the indefinite integral (with zero additive
%                    constant) of Chebyshev polynomials of the first kind
%                    from order 0 to order n evaluated at the points given
%                    in tau.
%
%INPUTS: tau     An NX1 or 1XN vector of values from tauStart to tauEnd, or
%                from -1 to 1 if tauStart and tauEnd are omitted, where one
%                wishes to evaluate the indefinite integrals
%                (antiderivatives) of the Chebyshev polynomials.
%        n       The non-negative integer maximum order of the Chebyshev
%                polynomials evaluated.
%tauStart,tauEnd The possible range of the inputs. If omitted, a  range of
%               -1 to 1 is assumed --the normal range for Chebyshev
%                polynomials. The option for mapping to a wider range is
%                useful when using Chebyshev polynomials for interpolation.
%                Note that the such polynomials are generally not useful
%                for interpolating much outside of the valid range.
%
%OUTPUTS: dT     An (n+1)XN matrix of the antiderivatives of the Chebyshev
%                polynomials from order 0 to n evaluated at each of the
%                values of tau. T(i,j) is the antiderivative of the (i-1)th
%                order Chebyshev polynomial evaluated at tau(j), taking
%                into account that the function has been mapped to the
%                range tauStart,tauEnd, if necessary.
%         T      Since the Chebyshev functions have to be computed to find
%                the indefinite integrals, an (n+2)XN matrix of the
%                Chebyshev functions can also be returned, if desired
%
%The recursion for finding the indefinite integral values can be obtained
%by examining how the function ChebyshevPolyIntCoeffs modifies the
%coefficients of a weighted series of Chebyshev polynomials when only a
%single 1 is present for each given order.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(tau);
tau=tau(:)';%Make tau into a column vector.
%Map the tau values to the -1 to 1 range, if necessary.
if(nargin>2)
    tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
end

%Get the Chebyshev polynomial function values.
T=ChebyshevPoly(tau,n+1);
intT=zeros(n+1,N);

%Deal with the n=0 case.
curN=0;
intT(curN+1,:)=T(curN+1+1,:);

if(n>0)
    %Deal with the n=1 case.
    curN=1;
    intT(curN+1,:)=0.5*(T(curN+1+1,:)/(curN+1)+0.5);
end

%Perform recusion for the rest of the values.
for curN=2:n
    intT(curN+1,:)=0.5*(T(curN+1+1,:)/(curN+1)-T(curN-1+1,:)/(curN-1));
end

%The constant to handle a mapping from -1->1 to some other range of
%parameters.
if(nargin>2)
    intT=intT*((tauEnd-tauStart)/2);
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
