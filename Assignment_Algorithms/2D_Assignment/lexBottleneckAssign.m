function [wAbs,wOrd,col4row]=lexBottleneckAssign(C)
%%LEXBOTTLENECKASSIGN Solve the lexicogrpahic bottleneck assignment problem
%          for a square cost matrix. A linear bottleneck assignment
%          problem can be formulated as
%          gain=minimize_phi max_i C_{i,phi(i)} where phi is a permutation
%                                               vector.
%          It can also be formulated as:
%          gain=minimize_{x} max_{i,j} C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{n}x_{i,j}=1 for all i
%          \sum_{i=1}^{n}x_{i,j}=1 for all j
%          x_{i,j}=0 or 1. Let the vector w be the vector of elements
%          C(i,phi(i)) sorted  in decreasing order. The lexicographic
%          lexicographic bottleneck problem solves the bottleneck problem
%          and produces a w that is lexicographically the minimum among all
%          possible w vectors. Note that gain=w(1).
%
%INPUTS: C An nXn cost matrix that does not contain any NaNs.
%
%OUTPUTS: wAbs The lexicographically minimal w vector. the elements are in
%              descending order.
%         wOrd The lexicographically minimal w vector given in terms of the
%              order statistics of the elements in C by magnitude. Thus, 0
%              is the minimum element in C, 1 is the second largest, etc.
%      col4row An nX1 vector where the entry in each element is an
%              assignment of the element in that row to a column. This is
%              equivalent to phi above.
%
%This function implements Algorithm 6.7 of Chapter 6.6 of [1]. A possibly
%faster technique is also described in the text of the Chapter.
%
%EXAMPLE 1:
%This is Example 6.2.4 in Chapter 6.6 of [1]. We compare the lexicographic
%bottleneck assignment problem solution to a generic bottleneck assignment
%problem solution.
% C=[8,9,3,2;
%    4,7,3,8;
%    0,8,10,4;
%    2,5,8,3];
% wAbsLex=lexBottleneckAssign(C);
% col4row=linBottleneckAssign(C);
% %Next, extract the  vector such that wAbs(k) is the kth largest cost in C
% %that is in the solution from linBottleneckAssign.
% n=size(C,1);
% wAbs=zeros(n,1);
% for i=1:n
%     j=col4row(i);
%     wAbs(i)=C(i,j);
% end
% wAbs=sort(wAbs,'descend');
% wAbsLex
% wAbs
%One sees that the ordered selected elements on wAbsLex have a lower
%lexicographic order than those of wAbs. wAbsLex is [5;3;2;0] and wAbs is
%[5;4;4;3].
%
%REFERENCES:
%[1] R. Burkard, M. Dell'Amico, and S. Martello,Assignment Problems.
%    Philadelphia: Society for Industrial and Applied Mathematics, 2009.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(C,1);

COrig=C;
[~,~,ic]=unique(C);
%C is now values from 0 to n^2-1, properly counting repeats at the same
%value.
C=reshape(ic-1,[n,n]);

%Solve the linear bottleneck assignment problem.
[col4row,~,z]=linBottleneckAssign(C);

%Allocate space for D
D=zeros(n,n);

F=false(n,n);
for i=1:n
    for j=1:n
        if(C(i,j)>z)
            F(i,j)=1;
        end
    end
end

for i=1:n
    for j=1:n
        if(F(i,j))
            D(i,j)=Inf; 
        else
            if(C(i,j)==z)
                D(i,j)=1;
            else
                D(i,j)=0;
            end
        end
    end
end

if(z>0)
    while(1)
        [col4row,~,~,u,v]=assign2D(D);

        z=z-1;
        if(z==0)
            break;
        end

        for i=1:n
            for j=1:n
                %The definitions of u and v used in assign2D are
                %opposite those used in Algorithm 6.7 of [1].
                if(D(i,j)>u(j)+v(i))
                    F(i,j)=1;
                end
            end
        end

        for i=1:n
            for j=1:n
                if(F(i,j))
                    D(i,j)=Inf; 
                else
                    if(C(i,j)==z)
                        D(i,j)=1;
                    else
                        D(i,j)=0;
                    end
                end
            end
        end
    end
end

wOrd=zeros(n,1);
wAbs=zeros(n,1);
for i=1:n
    j=col4row(i);
    wOrd(i)=C(i,j);
    wAbs(i)=COrig(i,j);
end

wAbs=sort(wAbs,'descend');
wOrd=sort(wOrd,'descend');
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
