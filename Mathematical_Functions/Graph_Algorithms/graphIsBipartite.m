function isBipartite=graphIsBipartite(A)
%%GRAPHISBIPARTITE Given a nondirectional graph represented by a symmetric
%                  adjacency matrix A, determine whether the graph is
%                  bipartite.
%
%INPUTS: A An nXn symmetric adjacency matrix. If element i is adjacent to
%          element j, then A(i,j) is nonzero. Otherwise, it is zero.
%          Diagonal elements are typically zero (if not, then the graph is
%          not bipartite).
%
%OUTPUTS: isBipartite This is true if the graph is bipartite and false
%                     otherwise.
%
%This function implements Algorithm B in Chapter 7 of [1], modified for an
%adjacency matrix. When very large, space graphs are used, then the
%adjacency matrix representation is not ideal.
%
%EXAMPLE:
%We consider four graphs. Two are bipartite and two are not.
%A1 and A2 are:
%A1:     A2:
%1--2    1--2
% \/      \/
% /\      /\
%3--4    3--6
%         \/
%         5
%and graphs A3 and A4 are:
%A3:        A4:
%1-----2    1-----2
%|\   /|    |\   /|
%| 5-6 |    | 5-6 |
%| | | |    | |\| |
%| 7-8 |    | 7-8 |
%|/   \|    |/   \|
%3-----4    3-----4
% A1=[0,1,0,1;
%     1,0,1,0;
%     0,1,0,1;
%     1,0,1,0];
% A2=[0,1,0,1,0;
%     1,0,1,0,0;
%     0,1,0,1,1;
%     1,0,1,0,1;
%     0,0,1,1,0];
% A3=[0,1,1,0,1,0,0,0;
%     1,0,0,1,0,0,0,0;
%     1,0,0,1,0,0,0,0;
%     0,1,1,0,0,0,0,0;
%     1,0,0,0,0,1,1,0;
%     0,0,0,0,1,0,0,1;
%     0,0,0,0,1,0,0,1;
%     0,0,0,0,0,1,1,0];
% A4=[0,1,1,0,1,0,0,0;
%     1,0,0,1,0,0,0,0;
%     1,0,0,1,0,0,0,0;
%     0,1,1,0,0,0,0,0;
%     1,0,0,0,0,1,1,1;
%     0,0,0,0,1,0,0,1;
%     0,0,0,0,1,0,0,1;
%     0,0,0,0,1,1,1,0];
% isBi1=graphIsBipartite(A1)
% isBi2=graphIsBipartite(A2)
% isBi3=graphIsBipartite(A3)
% isBi4=graphIsBipartite(A4)
%One will find that A1 and A3 are bipartite and A2 and A4 are not.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

%Allocate space.
link=zeros(n,1);

%Step B1
color=-1*ones(n,1);%In [1], it is indexed from 0, not 1; we index from 1.
w=n+1;%Only need to consider the ones originating an edge.

isBipartite=true;
while(w~=1)
    w=w-1;
    
    %Step B3
    if(color(w)>=0)
        continue;
    end
    color(w)=0;
    link(w)=0;
    s=w;
    
    while(1)
        %Step B4
        u=s;
        s=link(s);

        aList=find(A(u,:));

        for a=1:length(aList)
            v=aList(a);

            if(color(v)<0)
                color(v)=1-color(u);
                link(v)=s;
                s=v;
            elseif(color(v)==color(u))
                isBipartite=false;
                return;
            end
        end
        if(s==0)
            break;
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
