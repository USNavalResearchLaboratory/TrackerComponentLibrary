function intVal=betaIncomplete(x,a,b,normalize)
%%BETAINCOMPLETE Evaluate the incomplete beta function. This function
%  approximates the integral int_0^x t^a(1-t)^(b-1) dt if not normalzied.
%  If normalized, this integral is multiplied by 1/beta(a,b). Unlike
%  Matlab's built-in betainc, this implementation can handle -1<=x<=1.
%  That is, x can be negative, which as of Matlab R2022b is not allowed in
%  the betainc function. 
%
%INPUTS: x An real value -1<=x<=1.
%        a A positive scalar value.
%        b A positive scalar value.
%
%OUTPUTS: intVal The scalar value of the incomplete beta function.
%
%This function uses the identity between the incomplete beta function and
%the 2F1 hypergeometric function that is given in [1] and thus evaluates
%the incomplete beta function by calling hypergeometric2F1.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Incomplete Beta Function." From MathWorld--A
%    Wolfram Web Resource.
%    https://mathworld.wolfram.com/IncompleteBetaFunction.html
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(normalize))
    normalize=true;
end

if(a==0)
    intVal=1;
else
    intVal=(x^a/a)*hypergeometric2F1(a,1-b,a+1,x);
end

if(normalize)
    intVal=intVal/beta(a,b);
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
