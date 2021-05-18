function [theComp,recurData]=getNextLEQComposition(param1,q)
%%GETNEXTLEQCOMPOSITION Get the next composition of the integer n that has
%           q or fewer parts. A composition is a method of putting n
%           unlabeled items into a number of labeled slots. Compositions
%           are tuples of integers >=1 that sum to n.
%
%INPUTS: param If the first composition is desired, then param is n.
%              Otherwise, param is the output recurData from the previous
%              call to this function.
%            q The maximum number of parts. This value only needs to be
%              passed when generating the first composition q>1.
%
%OUTPUTS: theComp The tX1 composition or an empty matrix is the final
%                 composition has been reached.
%       recurData A data structure that can be passed back to this function
%                 to get subsequent compositions.
%
%This function implements Algorithm N in Problem 16 in Section 7.2.1.7 of
%[1].
%
%EXAMPLE:
%The example in the book.
% n=7;
% q=3;
% [theComp,recurData]=getNextLEQComposition(n,q);
% while(~isempty(theComp))
%     theComp'
%     [theComp,recurData]=getNextLEQComposition(recurData);
% end
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isstruct(param1))
    n=param1;
    q=min(q,n);
    a=zeros(q+1,1);%Allocate space
    
    %Step N1
    r=n;
    t=0;
    a(0+1)=0;
    recurData.q=q;
else
    recurData=param1;
    r=recurData.r;
    t=recurData.t;
    a=recurData.a;
    q=recurData.q;
end

if(isempty(a))
    theComp=[];
    return
end

%Step N2
while(r>=q)
    t=t+1;
    a(t+1)=q;
    r=r-q;
end

if(r>0)
    t=t+1;
    a(t+1)=r;
end

%Step N3
theComp=a(2:(t+1));

%Step N4
j=t;
while(a(j+1)==1)
    j=j-1;
end

if(j==0)
    recurData.a=[];
    return
end

%Step N5
a(j+1)=a(j+1)-1;
r=t-j+1;
t=j;

recurData.r=r;
recurData.t=t;
recurData.a=a;
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
