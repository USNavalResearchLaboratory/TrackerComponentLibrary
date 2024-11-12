function numVal=numFixedPointInvolutions(n)
%%NUMFIXEDPOINTINVOLUTIONS Determine the number of involutions of length n
%   where at least one element of the involution is in the same place it
%   started. An involution is a permution of items 1:n such that the
%   permutation is its own inverse.
%
%INPUTS: n The length of the involution n>=0. If 0, then 0 is returned.
%
%OUTPUTS: numVal The number of length n involutions that have at least one
%               fixed point.
%
%The sum implemented here is based on a principle of exclusion and
%inclusion. Suppose that out of n items, we choose one of the items to be
%fixed. That leads to  binomial(n,1)*numInvolutions(n-1) values. However,
%there are also duplicates of there being two or more things fixed there.
%After working out a principle of exclusion and inclusion to avoid double
%counting, we find that we sum (-1)^(k+1)*binomial(n,k)*numInvolutions(n-k)
%for k=1 to n where the flipping sign prevents double counting. Also, if n
%is odd, then this just returns numInvolutions(n) directly without summing,
%because it is not possible to have fixed-point-free involutions (pure
%involutions) when n is odd.
%
%EXAMPLE:
%In this example, we generate all of the involutions for a particular n and
%then count how many have fixed points. We show that the solution equals
%the value returned by this function.
% n=12;
% p=genAllInvolutions(n);
% pFixed=(1:n)';
% numInv=numInvolutions(n);
% sumVal=0;
% for k=1:numInv
%     sumVal=sumVal+any(pFixed==p(:,k));
% end
% sumVal
% numFixedPointInvolutions(n)
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(mod(n,2)==1)
    numVal=numInvolutions(n);
else
    numVal=0;
    for k=1:n
        numVal=numVal+(-1)^(k+1)*binomial(n,k)*numInvolutions(n-k);
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
