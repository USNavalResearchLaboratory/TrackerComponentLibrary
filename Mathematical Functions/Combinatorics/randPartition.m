function [thePartition,k,mult]=randPartition(n)
%%RANDPARTITION  Get a random partition of the integer n. A partition is a
%                set of other integers that sum to n. The partitions all
%                have equal probability.
%
%INPUTS: n A positive integer of which one desires a random partition.
%
%OUTPUTS: thePartition The random partition. This is a vector of integers
%                      that sum to n.
%                    k The number of unique digits in the partition.
%                 mult An nX1 vector where mult(i) is the number of times
%                      the digit i appears in the partition. Note that k
%                      and mult together uniquely specify the partition.
%
%The algorithm is based on the RANPAR function described in Chapter 10 of
%[1]. The algorithm has been modified to call the numberOfPartitions rather
%than using its own recursion to find the number of partitions from 1 to n.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Get the number of partitions for integers up through n.
p=numberOfPartitions(n,true);

%Get rid of the zero partition number.
p=p(2:end);

%Line 30 in the book
m=n;
k=0;
mult=zeros(n,1);

while(1)
    %Line 40 in the book
    z=rand(1)*m*p(m);
    d=0;

    quitLoop=false;
    while(1)
        if(quitLoop)
            break;
        end
        
        %Line 110 in the book
        d=d+1;
        i1=m;
        j=0;
        while(1)
            %Line 150 in the book.
            j=j+1;
            i1=i1-d;

            if(i1<0)
                break;
            elseif(i1==0)
                z=z-d;
                if(z<=0)
                    quitLoop=true;
                    break;
                else
                    break;
                end
            else
                z=z-d*p(i1);
                if(z<=0)
                    quitLoop=true;
                    break;
                else
                    continue;
                end
            end
        end
    end
    mult(d)=mult(d)+j;
    k=k+j;
    m=i1;
    
    if(m==0)
        break;
    end
end

%Build the partition including the repeats.
thePartition=zeros(k,1);

partIdx=1;
for idx=1:n
    numRep=mult(idx);
    thePartition(partIdx:(partIdx+numRep-1))=idx;
    partIdx=partIdx+numRep;
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
