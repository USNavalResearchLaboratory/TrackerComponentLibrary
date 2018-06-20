function isMonge=isMongeMatrix(A)
%%ISMONGEMATRIX Determine whether a matrix is a Monge matrix. An n1Xn2
%          matrix is a Monge matrix if A(i1,j1)+A(i2,j2)<=A(i1,j2)+A(i2,j1)
%          for all 1<=i1<i2<=n1 and 1<=j1<j2<=n2. Monge matrices have
%          special properties in some discrete optimization problems.
%           
%INPUTS: A An n1Xn2 matrix.
%
%OUTPUTS: isMonge This is true if A is a Monge matrix and false otherwise.
%
%EXAMPLE:
% A=[35, 25, 42, 88;
%    84, 63, 75, 96;
%    51,  9,  4, 25;
%    93, 16,  6, 18];
% isMonge=isMongeMatrix(A)
%A is a Monge matrix so isMonge will be 1.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n1=size(A,1);
n2=size(A,2);

isMonge=true;
for i2=2:n1
    for i1=1:(i2-1)
        for j2=2:n2
            for j1=1:(j2-1)
                if(A(i1,j1)+A(i2,j2)>A(i1,j2)+A(i2,j1))
                   isMonge=false;
                   return;
                end
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
