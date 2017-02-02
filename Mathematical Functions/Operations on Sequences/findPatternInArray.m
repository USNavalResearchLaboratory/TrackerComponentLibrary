function startIdx=findPatternInArray(pattern,array,maxTimes)
%%FINDPATTERNINARRAY Find all or a maximum number of occurrences of a
%                    pattern in an array of values.
%
%INPUTS: pattern A 1Xm or mX1 array holding elements of a pattern to be
%                found.
%        array   A 1Xn or nX1 array in which occurrences of the pattern are
%                to be found.
%       maxTimes An optional parameter specifying the maximum number of
%                occurrences of the pattern to find in array. If omitted or
%                an empty matrix is passed, then all occurrences of pattern
%                in array will be found.
%
%OUTPUTS: startIdx A numFoundX1 array of the atarting indices of the
%                occurrences of pattern in array, given in increasing
%                order. 
%
%The algorithm used is a linear time method given in [1]. It is the
%algorithm of Section 2 with the modification to the next array of Section
%4 to find all occurrences of the pattern. Not all of the suggestions for
%improving efficiency of Section 3 are used, because they would require
%lengthening the input arrays and also would require defining certain
%values the inputs could not take.
%
%EXAMPLE:
% array='baabbabbaabaabbaabbaabbabaa';
% pattern='aabba';
% startIdx=findPatternInArray(pattern,array)
%One finds the string at startIdx=[2;12;16;20]. Note that some of the
%solutions overlap.
%
%REFERENCES:
%[1] D. E. Knuth, J. H. Morris Jr., and V. R. Pratty, "Fast pattern
%    matching in strings," SIAM Journal on Computing, vol. 6, no. 2, pp.
%    323-350, Jun. 1977.
%
%August 2015, David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=length(pattern);
n=length(array);

if(m>n)
    %If the pattern is longer than the array, then a match is impossible.
    startIdx=[];
    return;
end

if(nargin<3||isempty(maxTimes))
   maxTimes=n-m+1;%The maximum
else
   maxTimes=min([n-m+1,maxTimes]);
end

%Allocate space for the maximum possible number of matches.
startIdx=zeros(maxTimes,1);
numFound=0;

%First, compute the elements of the next table. We are computing them for a
%length m+1 table where the last element does not match anything else so
%that we can get a resume index to use when finding multiple matches.
next=zeros(m+1,1);
j=1;
t=0;
next(1)=0;
while(j<m)
   while((t>0)&&pattern(j)~=pattern(t))
       t=next(t);
   end
   t=t+1;
   j=j+1;
   if(pattern(j)==pattern(t))
       next(j)=next(j);
   else
        next(j)=t;
   end
end

%Add in the value for j=m;
while((t>0)&&pattern(j)~=pattern(t))
   t=next(t);
end
t=t+1;
j=j+1;
next(j)=t;%This sets next(m+1);

%Next, run the actual algorithm using the next table.
j=1;
k=1;

%Consider the special case of j=1 --skip forward to the first possible
%match.
while(array(k)~=pattern(1))
    k=k+1;
    if(k>n)%Nothing matches the first character of the pattern.
        startIdx=[];
        return;
    end
end

while(1)
    while((j<=m)&&(k<=n))
        while((j>0)&&(array(k)~=pattern(j)))
           j=next(j);
        end
        k=k+1;
        j=j+1;
    end

    if(j<=m)%All of the matches have been found.
        startIdx=startIdx(1:numFound);
        return;
    end
    
    %With j>m, the leftmost match has been found in positions k-m through
    %k-1.
    numFound=numFound+1;
    startIdx(numFound)=k-m;
    j=next(m+1);
    
    if(numFound==maxTimes)
        startIdx=startIdx(1:numFound);
        return;
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
