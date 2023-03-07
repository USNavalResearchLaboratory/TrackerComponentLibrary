function intVal=monomialIntPyramid(alpha)
%%MONOMIALINTPYRAMID Evaluate the integral of prod_{i=1}^3 x(i)^alpha(i)
%  over the 3D pyramid with a square base with the peak vertex at (0,0,1)
%  and the base vertices at (1,-1,-1), (-1,-1,-1), (-1,1,-1), and (1,1,-1).
%
%INPUTS: alpha A 3X1 or 1X3 vector of the integer exponents of the
%              monomial term. All elements must be >=0.
%
%OUTPUTS: intVal The value of the specified integral.
%
%The solution was obtained via analytic integration and simplified to be in
%terms of beta functions. These types of explicit moment formulae are
%useful for testing cubature integration points.
%
%EXAMPLE:
%An analytic solution for alpha=[10;18;6] is 832273/941392221. Here, we
%show that this function produces the same result within resonable finite
%precision limits.
% alpha=[10;18;6];
% intSol=monomialIntPyramid(alpha);
% trueSol=832273/941392221;
% RelErr=(intSol-trueSol)/trueSol
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(all(mod(alpha(1:2),2)==0))
    a=alpha(1);
    b=alpha(2);
    c=alpha(3);

    intVal=(beta(1+c,3+a+b)-betaIncomplete(-1,1+c,3+a+b,false))/(2^(a+b)*(1+a)*(1+b));
else
    intVal=0;
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
