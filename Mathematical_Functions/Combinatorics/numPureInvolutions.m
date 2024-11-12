function numVal=numPureInvolutions(n)
%%NUMPUREINVOLUTIONS Determine the number of involutions of length n
%   where none of the elements of the involution is in the same place it
%   started. An involution is a permution of items 1:n such that the
%   permutation is its own inverse. Pure involutions (fixed-point free
%   involutions) are both involutions and derangements and they arise in
%   cryptanalysis, among other areas.
%
%INPUTS: n The length of the involution n>=0. If 0, then 0 is returned.
%
%OUTPUTS:numVal The number of length n involutions that have at least one
%               fixed point.
%
%This function just calls numInvolutions and subtracts
%numFixedPointInvolutions. If n is odd, then this just returns 0, since it
%is not possible to have pure involutions with n being odd.
%
%EXAMPLE:
%In this example, we generate all of the involutions for a particular n and
%then count how many have no fixed points. We show that the solution equals
%the value returned by this function.
% n=10;
% p=genAllInvolutions(n);
% pFixed=(1:n)';
% numInv=numInvolutions(n);
% sumVal=0;
% for k=1:numInv
%     sumVal=sumVal+(~any(pFixed==p(:,k)));
% end
% sumVal
% numPureInvolutions(n)
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(mod(n,2)==1)
    numVal=0;
else
    numVal=numInvolutions(n)-numFixedPointInvolutions(n);
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
