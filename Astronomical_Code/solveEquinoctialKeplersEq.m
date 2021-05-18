function F=solveEquinoctialKeplersEq(lambda,h,k)
%%SOLVEEQUINOCTIALKEPLERSEQ Solve the Equinoctial Kepler's Equation.
%               This is solving for F where lambda=F+h*cos(F)-k*sin(F);
%               This type of equation arises when using equinoctial orbital
%               elements.
%
%INPUTS: lambda The (real) left-hand side of the equation, typically the
%               mean longitude.
%          h,k  The real parameters multiplying the cosine and since terms,
%               typically the components of the eccentricity vector.
%
%OUTPUTS: F The real solution to lambda=F+h*cos(F)-k*sin(F).
%
%The algorithm is taken from Section 7.1 of [1]. The algorithm simply uses
%Newton's method to solve the equation. An arbitrarily upper limit of 50
%iterations has been set.
%
%REFERENCES:
%[1] D. A. Danielson, C. P. Sagovac, B. Neta, and L. W. Early,
%    "Semianalytic satellite theory," Mathematics Department, Naval
%    Postgraduate School, Monterey, CA, Tech. Rep., 1995. [Online].
%    Available: http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix= html&identifier=ADA531136
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

maxIter=50;

F=lambda;
for curIter=1:maxIter
    FOld=F;
    sinF=sin(F);
    cosF=cos(F);
    num=F+h*cosF-k*sinF-lambda;
    denom=1-h*sinF-k*cosF;
    F=F-num/denom;
    
    %If it converged to finite precision limits.
    if(abs(F-FOld)<=eps(F))
        break;
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
