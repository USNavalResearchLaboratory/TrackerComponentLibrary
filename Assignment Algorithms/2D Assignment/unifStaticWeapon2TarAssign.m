function x=unifStaticWeapon2TarAssign(V,q,M)
%%UNIFSTATICWEAP2TARASSIGN  Solve the uniform static weapon-to-target
%            assignment problem. This problem is
%            minimize Sum_{j=1}^NV(j)*q(j)^x(j)
%            subject to Sum_{j=1}^Nx(j)=M
%                       All x are non-negative integers.
%            This type of optimization cost function arises when one has M
%            weapons, all of which have an equal probability of destroying
%            target j, 1-q(j), and a positive real number V_j indicates
%            some type of preference between the targets. The optimization
%            determines how many weapons should be assigned to each target,
%            where x(j) is the number of weapons on target j.
%
%INPUTS: V An NX1 vector of multiplicative values, as described above,
%          where V(j) is the value for the jth target. These should be
%          non-negative.
%        q An NX1 vector of values that are raised to the x values to be
%          found. These should be between 0 and 1.
%        M The required total sum of the x values found in the
%          optimization.
%
%OUTPUTS: x An NX1 vector indicating how many of the M weapons should be
%           assigned to the targets. This is the optimal solution.
%
%The algorithm is from [1], the global optimality of the algorithm is
%proven. The algorithm is also described in Chapter 3 of [2]. The algorithm
%is also known as the maximal marginal return (MMR) algorithm.
%
%REFERENCES:
%[1] G. G. den Broeder Jr., R. E. Ellison, and L. Emerling, "On optimum
%    target assignments," Operations Research, vol. 7, no. 3, pp. 322-326,
%    May - Jun. 1959.
%[2] P. M. Pardalos and L. Pitsoulis, Nonlinear Assignment Problems.
%    Dordrecht: Kluwer Academic Publishers, 2000.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(V);

%Allocate space for the return variables.
x=zeros(N,1);

%Allocate space
S=V;
invQ=1-q;

for weaponIndex=1:M
    %Find the target to which the weapon will do the most damage.
    [~,k]=max(S.*invQ);
    %Assign the weapon to the target.
    x(k)=x(k)+1;
    S(k)=S(k)*q(k);
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
