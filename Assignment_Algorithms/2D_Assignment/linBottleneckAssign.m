function [col4row,row4col,gain]=linBottleneckAssign(C,solSel)
%%LINBOTTLENECKASSIGN Solve the linear bottleneck assignment problem for a
%          square cost matrix. The problem being solved can be formulated
%          as
%          gain=minimize_phi max_i C_{i,phi(i)} where phi is a permutation
%                                               vector.
%          It can also be formulated as:
%          gain=minimize_{x} max_{i,j} C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{n}x_{i,j}=1 for all i
%          \sum_{i=1}^{n}x_{i,j}=1 for all j
%          x_{i,j}=0 or 1.
%
%INPUTS: C An nXn cost matrix that does not contain any NaNs.
%   solSel Multiple optimal solutions typically exist to this problem. This
%          parameter specifies which solution to choose. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Choose
%            a feasible solution without any particular property.
%          1 Choose the optimal solution that also minimizes
%            sum_i C_{i,phi(i)}.
%          2 Choose the optimal solution that also maximizes
%            sum_i C_{i,phi(i)}.
%
%OUTPUTS: col4row An nX1 vector where the entry in each element is an
%                 assignment of the element in that row to a column. This
%                 is equivalent to phi above.
%         row4col A nX1 vector where the entry in each element is an
%                 assignment of the element in that column to a row.
%            gain The value of the cost function of the linear bottleneck
%                 assignment problem at the optimal point.
%
%This function implements Algorithm 6.1 of Chapter 6.2 of [1]. This is an
%algorithm based on thresholding. The choice of which solution is obtained
%in the end depends on whether one uses maxCardBipartMatching in the final
%step, or whether one modifies the cost matrix in the final step according
%to G (forbidden assignments become either Inf or -Inf) and then uses
%assign2D to solve the problem.
%
%EXAMPLE:
%This is example 6.5 from Chapter 6.2 of [1].
% C=[8,2,3,3;
%    2,7,5,8;
%    0,9,8,4;
%    2,5,6,3];
% [col4row,row4col,gain]=linBottleneckAssign(C)
% %Verify the value of the cost function for the permutation.
% n=size(C,1);
% costVal=-Inf;
% for curRow=1:n
%     cij=C(curRow,col4row(curRow));
%     if(cij>costVal)
%         costVal=cij;
%     end
% end
% costVal
%One will see that gain=5 and costVal agrees with the gain. The cost value
%is five even though the col4row assignment is not identical to the in the
%book. Solution variations depend on variations to the algiruthm used for
%the perfect matching. Here, maxCardBipartMatching is used.
%
%REFERENCES:
%[1] R. Burkard, M. Dell'Amico, and S. Martello,Assignment Problems.
%    Philadelphia: Society for Industrial and Applied Mathematics, 2009.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(solSel))
    solSel=0;
end

n=size(C,1);

c0=min(C(:));
c1=max(C(:));

if(c0==c1)
    col4row=(1:n).';
    row4col=(1:n).';
    gain=c0; 
    return;
end

sel=C<c1&C>c0;
CSel=C(sel);
if(isempty(CSel))
    %This instance can occur, for example, if the matrix contains only two
    %values, such as if C is the identity matrix.
    cStar=c1;
else
    while(~isempty(CSel))
        cStar=medianHi(CSel(:));
        [c0,c1]=feasCheck(C,cStar,c0,c1);

        sel=C<c1&C>c0;
        CSel=C(sel);
    end
end

if(cStar~=c0)
    [~,c1]=feasCheck(C,c0,c0,c1);
end

G=(C<=c1);
switch(solSel)
    case 0
        [col4row,row4col]=maxCardBipartMatching(G);
    case 1
        CNew=Inf(n,n);
        CNew(G)=C(G);
        [col4row,row4col]=assign2D(CNew,false);
        if(isempty(col4row))%If the cost wouldn't have been finite.
            [col4row,row4col]=maxCardBipartMatching(G);
        end
    case 2
        CNew=-Inf(n,n);
        CNew(G)=C(G);
        [col4row,row4col]=assign2D(CNew,true);
        if(isempty(col4row))%If the cost wouldn't have been finite.
            [col4row,row4col]=maxCardBipartMatching(G);
        end
    otherwise
        error('Unknown solution selected.')
end
gain=c1;

end

function [c0,c1]=feasCheck(C,cStar,c0,c1)

G=(C<=cStar);
col4row=maxCardBipartMatching(G);

if(all(col4row~=0))
    %It is feasible.
    c1=cStar;
else
    %It is infeasible.
    c0=cStar;
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
