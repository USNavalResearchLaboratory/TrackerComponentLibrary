function val=solveJosephusProblem(n,k)
%%SOLVEJOSEPHUSPERMUTATION This solves the Josephus problem. Given a circle
% of n items, remove every kth item (going around the circle, so after
% removing 1 item, there are n-1 items) until there is only 1 item left.
% The value is the index of the remaining item. It coincides with the value
% in the final element of a Josephus permutation.
%
%INPUTS: n The total number of items.
%        k The kth item is taken each time. This is that number. Note that
%          k can be bigger than n, because it loops back.
%
%OUTPUTS: val The integer solution fo the Josephus problem. This is a value
%             from 1 to n.
%
%The origin and properties of the Josephus problem are discussed in detail
%in [1] from which one can come up with a solution method using addition
%and modulo operations. The problem is discusse din general in [2].
%
%EXAMPLE:
%These two are from the example Josephus permtuations on page 15 of [1].
%The solutions are 12 and 4:
% val=solveJosephusProblem(12,16805)
% val=solveJosephusProblem(14,167565)
%
%REFERENCES
%[1] P. Schumer, "The Josephus problem: Once more around," Mathematics
%    Magazine, vol. 75, no. 1, pp. 12-17, Feb. 2002.
%[2] Weisstein, Eric W. "Josephus Problem." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/JosephusProblem.html
%
%February 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=0;
for i=1:n
    val=mod(val+k,i);
end

val=val+1;

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
