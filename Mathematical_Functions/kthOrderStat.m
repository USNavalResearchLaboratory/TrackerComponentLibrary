function [retVal,A]=kthOrderStat(A,k,dim)
%%KTHORDERSTAT Given an unsorted real vector a, return the kth smallest
%              element of the vector. This is done without sorting the
%              vector, though the output vector a will be somewhat more
%              sorted compared to the input.
%
%INPUTS: A A real vector or matrix of values, whereby the kth order
%          statistic is taken over dimension dim. 
%        k The selected order of the value to return. 1<=k<=n, where n is
%          the size of dimension dim. k=1 means the smallest element, k=n
%          means the largest element.
%      dim An optional parameter specifying the dimensions across which
%          the kth order statistics are found. If this value is omitted or
%          an empty matrix is passed, then the order statistics are
%          computed across the first non-singleton dimension of x.
%
%OUTPUTS: retVal The kth smallest element of A. This has the same shape as
%                A with dimension dim reduced to a singleton. If an empty
%                matrix is passed for A, then retVal is an empty matrix.
%              A The modified form of A. When performing multiple calls
%                with different values of k, it will typically be faster if
%                the modified A is passed back.
%
%This function implements the Floyd-Rivest algorithm, discussed
%theoretically in [1] and explicitely in [2]. However, the suggested
%adjustments involving exponentiation and square roots when dealing with
%very large vectors has been omitted. This is O(n) complexity, whereas
%sorting can be O(n*log(n)).
%
%EXAMPLE:
%Here, we demonstrate that the result is the same as
% numRuns=1000;
% numEls=20;
% for curRun=1:numRuns
%     a=randn(numEls,1);
%     k=randi(numEls);
% 
%     kthVal=kthOrderStat(a,k);
% 
%     a=sort(a,'ascend');
%     valAlt=a(k);
%     %Test that it is equivalent to sorting and choosing the kth element
%     assert(valAlt==kthVal)
% end
%There should be no error (due to the assertion) because the two techniques
%are equivalent.
%
%REFERENCES:
%[1] R. W. Floyd and R. L. Rivest, "Expected time bounds for selection,"
%    Communications of the ACM vol. 18, no. 3, pp. 165-172, Mar. 1975.
%[2] R. W. Floyd and R. L. Rivest, "Algorithm 489: The algorithm SELECT for
%    finding the ith smallest of n elements [M1]," Communications of the
%    AES, vol. 18, no. 3, p. 173, Mar. 1975.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If an empty matrix is passed.
if(isempty(A))
    retVal=[];
    return;
end

aSize=size(A);
nDims=length(aSize);

%Select the first non-singleton dimension if dim is not provided.
if(nargin<3||isempty(dim))
    dim=find(aSize>1,1);
    if(isempty(dim))
        dim=1;
    end
end

IdxBeforeDim=prod(aSize(1:(dim-1)));
IdxAfterDim=prod(aSize((dim+1):nDims));
n=aSize(dim);

if(k>n)
    error('k cannot be >n');
end

A=reshape(A,[IdxBeforeDim,n,IdxAfterDim]);

retVal=zeros(IdxBeforeDim,1,IdxAfterDim);
for curBefore=1:IdxBeforeDim
    for curAfter=1:IdxAfterDim
        a=A(curBefore,:,curAfter);

        l=1;
        r=n;
        while(r>l)
            t=a(k);

            i=l;
            j=r;

            %Swap a(l) and a(k)
            temp=a(k);
            a(k)=a(l);
            a(l)=temp;

            if(a(r)>t)
                %Swap a(r) and a(l)
                temp=a(r);
                a(r)=a(l);
                a(l)=temp;
            end

            while(i<j)
                %Swap a(i) and a(j)
                temp=a(i);
                a(i)=a(j);
                a(j)=temp;

                i=i+1;
                j=j-1;
                while(a(i)<t)
                    i=i+1;
                end

                while(a(j)>t)
                    j=j-1;
                end
            end

            if(a(l)==t)
                %Swap a(l) and a(j)
                temp=a(l);
                a(l)=a(j);
                a(j)=temp;
            else
                j=j+1;
                %Swap a(j) and a(r)
                temp=a(j);
                a(j)=a(r);
                a(r)=temp;
            end

            if(j<=k)
                l=j+1;
            end

            if(k<=j)
                r=j-1;
            end
        end

        retVal(curBefore,1,curAfter)=a(k);
        A(curBefore,:,curAfter)=a;
    end
end

A=reshape(A,aSize);

aSize(dim)=1;
retVal=reshape(retVal,aSize);

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
