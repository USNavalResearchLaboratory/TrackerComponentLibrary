function An=menageNumber(n)
%%MENAGENUMBER Get the nth ménage number. This is the number of ways that n
%       married couples can be seated around a circular table such that
%       there is always one man between two women and no man is next to his
%       own wife.
%
%INPUTS: n The integer number of couples that are to be seated. n>=2.
%
%OUTPUTS: An The number of seating arrangements such that each man is
%            always seated between two women.
%
%The problem is solved using Laisant's recurrence formula, which is derived
%in Chapter 8 of [1].
%
%REFERENCES:
%[1] H. Dörrie, 100 Great Problems of Elementary Mathematics: Their History
%    and Solution. New York: Dover Publications, Inc., 1965.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n<2)
    error('n must be >=2.');
end

%Special cases
if(n==2)
    An=0;
    return;
elseif(n==3)
    An=1;
    return;
end

Ak2=0;%A2
Ak1=1;%A3
for k=4:n
    ANext=(-4*(-1)^k+Ak2*k+Ak1*k*(k-2))/(k-2);
    Ak2=Ak1;
    Ak1=ANext;
end
An=Ak1;

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
