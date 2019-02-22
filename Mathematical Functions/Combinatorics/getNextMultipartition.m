function [vOut,data]=getNextMultipartition(data)
%%GETNEXTMULTIPARTITION Get the next multipartition of a set of elements
%           with possible repeats, or get the first multipartition in the
%           sequence. The multipartitions are given in decreasing
%           lexicographic order.
%
%INPUTS: data If the first multipartition is desired, then this is an mX1
%             or 1Xm vector containing how many times each element in the
%             set to be partitioned is repeated. For example, if one wishes
%             to partition the set (1,2,3,4), then this is [1;1;1;1]. For
%             the set [1;2;2;2;3;3;4], then this is [1;3;2;1]. Otherwise,
%             if one wishes to get the next multipartition in the sequence,
%             then this is the data variable returned from the previous
%             call to this function.
%
%OUTPUTS: vOut This is the current multipartition, or an empty vector if
%              one has passed the final multipartition in the sequence. The
%              multipartition is an mXnumParts matrix. The number of
%              columns is the number of parts of the partition. Each column
%              has the number of repeats of each item (corresponding to the
%              positions in the original data passed). See the examples
%              below for more explanation.
%         data This is a structure that can be passed to this function to
%              get the next multipartition in the sequence.
% 
%This function implements Algorithm M of Chapter 7.2.1.5 of [1].
%
%As an example, there are nine multipartitions of the set (1,1,2,2). These
%are
%(1,1,2,2);  (1,1,2),(2);  (1,1),(2,2);  (1,1),(2),(2); (1,2,2),(1);
%(1,2),(1,2);  (1,2),(1),(2);  (1),(1),(2,2); (1),(1),(2),(2)
%However, represented in terms of the output vOut, we have 9 matrices:
% [2;  [2,0; [2,0; [2,0,0;  [1,1;  [1,1;  [1,1,0;  [1,1,0;  [1,1,0,0
%  2]   1,1]  0,2]  0,1,1]   2,0]   1,1]   1,0,1]   0,0,2]   0,0,1,1]
%The example below shows how to use this function and how to display the
%outputs in a more recognizable form as partitions.
%
%EXAMPLE:
%Here, we display all multipartitions of 1123, which means that n=[2;1;1].
% vTotal=[1;1;2;3];
% vUnique=unique(vTotal);
% n=[2;1;1];
% m=length(n);
% [vOut,data]=getNextMultipartition(n);
% while(~isempty(vOut))
%     %Display the current multipartition.
%     %We will build a string showing the partitions. of [1,1,2,3].
%     strDisp=[];
%     
%     numParts=size(vOut,2);
%     for curPart=1:numParts
%         strDisp=[strDisp,sprintf('(')];
%         for k=1:m
%             for i=1:vOut(k,curPart)
%                 strDisp=[strDisp,sprintf('%i',vUnique(k))];
%             end
%         end
%         strDisp=[strDisp,sprintf(')')];
%     end
%     disp(strDisp)
%     
%     %Get the next one.
%     [vOut,data]=getNextMultipartition(data);
% end
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If the first multipartition is desired.
if(~isstruct(data))
    n=data;
    data=[];

    m=length(n);
    nSum=sum(n);
    
    data.nSum=nSum;
    data.m=m;

    %Allocate space
    f=zeros(nSum+1,1);
    c=zeros(m*nSum+1,1);
    u=zeros(m*nSum+1,1);
    v=zeros(m*nSum+1,1);

    %Step M1
    c(1:m)=1:m;
    u(1:m)=n(1:m);
    v(1:m)=n(1:m);
    f(0+1)=0;
    a=0;
    l=0;

    f(1+1)=m;
    b=m;
else
    %If the algorithm is finished.
    if(data.algFinished)
       vOut=[];
       return;
    end
    
    f=data.f;
    c=data.c;
    u=data.u;
    v=data.v;
    a=data.a;
    b=data.b;
    l=data.l;
    m=data.m;
    nSum=data.nSum;
end

while(1)
    %Step M2
    j=a;
    k=b;
    x=0;

    while(j<b)
        u(k+1)=u(j+1)-v(j+1);

        if(u(k+1)==0)
            x=1;
            j=j+1;
        elseif(x==0)
            c(k+1)=c(j+1);
            v(k+1)=min(v(j+1),u(k+1));
            x=u(k+1)<v(j+1);
            k=k+1;
            j=j+1;
        else
            c(k+1)=c(j+1);
            v(k+1)=u(k+1);
            k=k+1;
            j=j+1;
        end
    end

    %Step M3
    if(k>b)
        a=b;
        b=k;
        l=l+1;
        f(l+1+1)=b;
        continue;%Return to M2
    end

    %Step M4, visit the partition
    numOut=0;
    vOut=zeros(m,nSum);%Allocate the maximum possible space needed.
    for k=0:l
        numOut=numOut+1;
        for j=f(k+1):(f(k+1+1)-1)
            vOut(c(j+1),numOut)=v(j+1);
        end
    end
    vOut=vOut(:,1:numOut);
    
    %Step M5
    while(1)
    
        j=b-1;
        while(v(j+1)==0)
            j=j-1; 
        end

        if(~(j==a&&v(j+1)==1))
            v(j+1)=v(j+1)-1;

            for k=(j+1):(b-1)
               v(k+1)=u(k+1); 
            end

            %Return to M2 to visit the next one, so return from this
            %function.
            data.f=f;
            data.c=c;
            data.u=u;
            data.v=v;
            data.a=a;
            data.b=b;
            data.l=l;
            data.algFinished=false;
            return
        end

        %Step M6
        if(l==0)%All multiset partitions have been found.
            %There is no need to update the other things in data, because
            %there are no more multipartitions to find.
            data.algFinished=true;
            return;
        end

        l=l-1;
        b=a;
        a=f(l+1);
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
