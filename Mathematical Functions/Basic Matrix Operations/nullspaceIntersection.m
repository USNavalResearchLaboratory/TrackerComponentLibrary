function Y=nullspaceIntersection(A,B,algorithm)
%%NULLSPACEINTERSECTION Given an mXn matrix A and a pXn matrix B, find a
%               basis for the intersection of the nullspaces of the
%               matrices. That is, the output Y is a vectors whose columns
%               are such that A*Y(:,i)=0 and B*Y(:,i) =0 for all i columns.
%               If there is no nullspace intersection, then return an empty
%               matrix.
%
%INPUTS: A An mXn matrix.
%        B A pXn matrix.
% algorithm An optional parameter specifying which algorithm to use.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Algorithm 6.4.2 in Chapter 6.4.2 of [1].
%          1 In 6.4.2, it is noted that the common nullspace can be found
%            by finding the nullspace of [A;B]. This algorithm just passes
%            the augmented matrix to the null function.
%
%OUTPUTS: Y A basis for the common nullspace of A and B. If there is no
%           common nullspace, then an empty matrix is returned.
%
%In algorithm 0, when determining the rank  of the matrix C in the
%algorithm, matrix B is used in the scaling so as to handle the case where
%the rank should be 0 because B and one set of singular vectors of A are
%orthogonal.
%
%EXAMPLE:
%This just gives an example where both approaches give the same result,
%except the sign of the answer is different.
% A=[16,  2,  3,  13;
%     5, 11, 10,   8;
%     9,  7,  6,  12];
% B=[9,  7,  6, 12;
%    4, 14, 15,  1];
% nullspaceIntersection(A,B,0)
% nullspaceIntersection(A,B,1)
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0
        n=size(A,2);
        [~,SA,VA]=svd(A,0);

        r=sum(diag(SA)>eps()*norm(A,1));
        if(r<n)
            C=B*VA(:,r+1:n);
            [~,SC,VC]=svd(C,0);
            
            %We use B instead of C for the scaling in determining the rank,
            %to get the correct scaling for the case where the rank might
            %end up to be 0, because C is effectively an all zero matrix (B
            %and VA being orthogonal).
            q=sum(diag(SC)>eps()*norm(B,1));
            if(q<n-r)
                Y=VA(:,(r+1):n)*VC(:,(q+1):(n-r));
            else
                Y=[];%There is no common nullspace.
            end
        else
            %If A has no nullspace.
            Y=[]; 
        end
    case 1
        Y=null([A;B]);
    otherwise
        error('Unknwon algorithm specified.')
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
