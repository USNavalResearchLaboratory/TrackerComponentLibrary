function thePerms=genAllConstrainedPerms(n,t,preallocate)
%%GENALLCONSTRAINEDPERMS This function generates all permuations with a
%                user-specified constraint on prefixes. The permutation
%                must pass the tests t(x(1))==true, t(x(1:2))==true,
%                t(x(1:3))==true, etc. to be accepted. Thus, if any of the
%                tests on the prefixes fail, all permutations with those
%                prefixes are discarded. The permutations are visited in
%                lexicographic order.
%
%INPUTS: n The positive integer dimensionality of the permutations desired.
%        t The optional constraint function handle. If omitted, all
%          permutations are generated. This function should take partial
%          permutations as inputs and return 0 or 1.
% preallocate An optional boolean variable indivating whether the full
%          factorial(n) (maximum) output matrix should be preallocated for
%          the return variable. The default if this parameter is omitted or
%          an empty matrix is passed is true. A reason to not preallocate
%          would be if t greatly reduces the number of permutations that
%          need to be visited and factorial(n) is very large.
%
%OUTPUTS: thePerms The nXnumPerms set of permutations. The maximum number
%                  of permutations is factorial(n) (whicb is preallocated),
%                  but the constraints in t might case fewer to be
%                  returned.
%
%This function implements Algorithm X of Chapter 7.2.1.2 of [1].
%
%EXAMPLE:
%Here, we generate all permutations of 4 items:
% allPerms=genAllConstrainedPerms(4)
% %Now we generate all permutations of 4 items such that the sum of the
% %first two entries is less than or equal to four.
% t=@(x)(sum(x(1:min(length(x),2)))<=4);
% permsConst=genAllConstrainedPerms(4,t)
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(t))
    t=@(x)(1);
end

if(nargin<3||isempty(preallocate))
    preallocate=true;
end

numPerms=factorial(n);%The maximum possible number of permutations.
if(preallocate==true)
    thePerms=zeros(n,numPerms);
else
    thePerms=zeros(n,1);
end

a=zeros(n,1);
u=zeros(n,1);

curPerm=1;

%Step X1
l=[1:n,0];
k=1;

while(1)
    %Step X2, enter level k
    p=0;
    q=l(0+1);

    while(1)
        %Step X3
        a(k)=q;
        testVal=t(a(1:k));
        if(testVal==true&&k==n)
            thePerms(:,curPerm)=a;
            curPerm=curPerm+1;
            %Goto X6
        elseif(testVal==true&&k~=n)
            %Step X4
             u(k)=p;
             l(p+1)=l(q+1);
             k=k+1;
             break;%Go to X2
        elseif(testVal==false)
            %Step X5
            p=q;
            q=l(p+1);
            if(q~=0)
                continue;%Go to X3
            end
        end

        while(1)
            %Step X6
            k=k-1;
            if(k==0)
                %Shrink to fit
                thePerms=thePerms(:,1:(curPerm-1));
                %return
                return;
            end

            p=u(k);
            q=a(k);
            l(p+1)=q;
        
            %Step X5
            p=q;
            q=l(p+1);
            if(q~=0)
                break;%Go to X3
            end
        end
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
