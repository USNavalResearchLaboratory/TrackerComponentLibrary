function [pPart,cPart,recurData]=getNextBinPartition(param1)
%%GETNEXTBINPARTITION Get the next partition of an integer n into binary
%           parts. That is, n=sum(pPart.*cPart), where pPart is a vector of
%           powers of 2 and cPart are integer coefficients >=1.
%
%INPUTS: param1 To get the first partition in the sequence, run this
%               function with param1=n, the integer>=1 to be partitioned.
%               To get subsequent values, use param1=recurData.
%
%OUTPUTS: pPart, cPart The numPartX1 vectors expressing the current
%               partition with n=sum(pPart.*cPart), pPart being powers of 2
%               and cPart being integers>=1. The number of parts will vary.
%               When the final binary partition has been passed, empty
%               matrices will be returned.
%     recurData A data structure that can be passed back to this function
%               to get subsequent binary partitions.
%
%This function implements Problem 62 of Section 7.2.1.4 of [1].
%
%EXAMPLE:
%Here, we just verify that the partitioned values actually equal the
%desired value.
% n=51;
% [pPart,cPart,recurData]=getNextBinPartition(n);
% pLenOrig=length(recurData.p);
% while(~isempty(recurData))
%     assert(sum(pPart.*cPart)==n)
%     assert(length(recurData.p)==pLenOrig)
%     [pPart,cPart,recurData]=getNextBinPartition(recurData);
% end
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If this is the first partition.
if(~isstruct(param1))
    n=param1;
    maxParts=floor(log2(n))+1;
    %Indexation of c and p in the book is from 0.
    c=zeros(maxParts,1);
    p=zeros(maxParts,1);

    %Step K1
    c(0+1)=0;
    p(0+1)=0;
    p(1+1)=1;
    
    if(mod(n,2)==0)%If n is even
        c(1+1)=n;
        t=1;
    else
        k=1;
        q=(n-1)/2;
        while(mod(q,2)==0)
            q=q/2;
            k=k+1;
        end
        %Now (n-1)=q*2^k with q being odd.
        c(1+1)=1;
        c(2+1)=q;
        p(2+1)=2^k;
        t=2;
    end
    
    %Step even visit.
    pPart=p(2:(t+1));
    cPart=c(2:(t+1));
    
    recurData.p=p;
    recurData.c=c;
    recurData.t=t;
    recurData.evenVisit=true;
    return
else
    recurData=param1;
    
    p=recurData.p;
    c=recurData.c;
    t=recurData.t;
    evenVisit=recurData.evenVisit;
end

if(evenVisit)
    %Step K3, change the largest part.
    if(c(t+1)==1)
        %Split the largest part.
        if(p(t+1)~=2*p(t-1+1))
            c(t+1)=2;
            p(t+1)=p(t+1)/2;
        else
           c(t-1+1)=c(t-1+1)+2;
           t=t-1;
        end
    else%c(t)>1
        %Merge the two largest parts.
        if(c(t+1)==2)
            c(t+1)=1;
            p(t+1)=2*p(t+1);
        else
            c(t+1)=c(t+1)-2;
            c(t+1+1)=1;
            p(t+1+1)=2*p(t+1);
            t=t+1;
        end
    end

    %%Step odd visit.
    pPart=p(2:(t+1));
    cPart=c(2:(t+1));
    recurData.p=p;
    recurData.c=c;
    recurData.t=t;
    recurData.evenVisit=false;
    return
end

%Step K5, change the next largest part. 9 cases.
if(mod(c(t+1),2)==1)%If c(t) is odd.
    if(t==1)%1a
        pPart=[];
        cPart=[];
        recurData=[];
        return
    elseif(c(t-1+1)==1)
        if(p(t-1+1)==2*p(t-2+1))%1b1
            c(t-2+1)=c(t-2+1)+2;
            c(t-1+1)=c(t+1);
            p(t-1+1)=p(t+1);
            t=t-1;
        else%1b2
            c(t-1+1)=2;
            p(t-1+1)=p(t-1+1)/2;
        end
    elseif(c(t-1+1)==2)
        if(p(t+1)==2*p(t-1+1))%1c1
            c(t-1+1)=c(t+1)+1;
            p(t-1+1)=p(t+1);
            t=t-1;
        else%1c2
            c(t-1+1)=1;
            p(t-1+1)=2*p(t-1+1);
        end
    else%c(t-1)>2
        if(p(t+1)==2*p(t-1+1))%1d1
            c(t-1+1)=c(t-1+1)-2;
            c(t+1)=c(t+1)+1;
        else%1d2
            c(t+1+1)=c(t+1);
            p(t+1+1)=p(t+1);
            c(t+1)=1;
            p(t+1)=2*p(t-1+1);
            c(t-1+1)=c(t-1+1)-2;
            t=t+1;            
        end
    end
else%If c(t) is even.
    if(p(t+1)==2*p(t-1+1))%2a
        c(t+1)=c(t+1)-1;
        c(t-1+1)=c(t-1+1)+2;
    else%2b
        c(t+1+1)=c(t+1)-1;
        p(t+1+1)=p(t+1);
        c(t+1)=2;
        p(t+1)=p(t+1)/2;
        t=t+1;
    end
end

%Step even visit.
pPart=p(2:(t+1));
cPart=c(2:(t+1));

recurData.p=p;
recurData.c=c;
recurData.t=t;
recurData.evenVisit=true;

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
