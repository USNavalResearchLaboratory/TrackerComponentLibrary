function numPart=numberOfPartitions(n,returnAll)
%%NUMBEROFPARTITIONS The number of ways of partitioning the positive
%                    integer n. This is the number of ways of writing n as
%                    a sum of other positive integers. This is sometimes
%                    called partition function p. If one only wants the
%                    number of possible ways of partitioning n into m
%                    parts, then the function numMPartitions should be
%                    used.
%
%INPUTS: n The positive integer, n>0, that is to be considered for
%          partitioning.
% returnAll An optional parameter specifying what should be returned.
%          Possible values are
%          false (The default if omitted or an empty matrix is passed) Only
%                 the number of partitions of n is returned.
%           true  A vector of all of the numbers of partitions from n=0 to
%                 n are returned. The number of partitions of a positive
%                 integer 0<=i<=n is thus numPart(i+1);
%
%OUTPUTS: numPart The number of possible partitions of the integer n, if
%                 returnAll is omitted or is false, or an (n+1)X1 vector of
%                 the number of partitions of all integers from 0 to n.
%
%For n<=250, the algorithm used is the recursion taken from Chapter
%7.2.1.4 of [1]. The recursion is also mentioned in [2]. However, as n gets
%larger, the recursion suffers from a loss of precision that can lead to
%significant inaccuracy at higher values of n. For example, at n=2000, the
%result is about 265 times too large. Thus, for n>250, an algorithm based
%on that given in the first half of the RANPAR algorithm of Chapter 10 of
%[3] is used. That algorithm is slower, when implemented in Matlab, but is
%much more robust to finite precision errors at large values.
%
%The amount of memory used for the computation scales linearly with n,
%because all past values of p have to be found during the recursion. For
%very large values of n, finite precision errors will eventually become
%large.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%[2] Weisstein, Eric W. "Partition Function P." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/PartitionFunctionP.html
%[3] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(returnAll))
    returnAll=false;
end

if(n<=250)
    %This is the algorthm of Knuth.
    p=zeros(n+1,1);
    p(0+1)=1;

    for curN=1:n
        k=1;
        signVal=1;
        while(1)
            %These are generalized pentagonal numbers.
            k1=(3*k^2-k)/2;
            if(curN-k1<0)
                break;
            end

            p(curN+1)=p(curN+1)+signVal*p(curN-k1+1);

            k2=(3*k^2+k)/2;
            if(curN-k2<0)
                break;
            end

            p(curN+1)=p(curN+1)+signVal*p(curN-k2+1);

            k=k+1;
            signVal=-signVal;
        end
    end
    
    if(returnAll)
        numPart=p;
    else
        numPart=p(end);
    end
else
    %This is the algorithm of Nijenhuis and Wilf.
    p=zeros(n,1);
    p(1)=1;
    m=1;
    if(n~=1)
        for i=m:n
            iSum=0;
            for d=1:i
                is=0;
                i1=i;

                while(1)
                    i1=i1-d;

                    if(i1<0)
                        break;
                    elseif(i1==0)
                        is=is+1;
                        break;
                    else
                        is=is+p(i1);
                    end
                end

                %Line 22
                iSum=iSum+is*d;
            end

            %Line 21
            p(i)=iSum/i;
        end
    end
    
    if(returnAll)
        numPart=[1;p];
    else
        numPart=p(end);
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
