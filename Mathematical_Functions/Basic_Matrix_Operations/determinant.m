function val=determinant(A,algorithm)
%%DETERMINANT This computes the determinant of a real matrix without using
%             any built-in Matlab matrix decompositions. The determinant is
%             computed using an LU decomposition or for matrices 5D and
%             below, if no algorithm is specified, then an explicit formula
%             is used. Alternatively, if the matrix is known to be
%             symmetric and positive definite, an algorithm using a
%             Cholesky decomposition can be used.
%
%INPUTS: A An nXn matrix.
% algorithm An optional parameter specifying the algorithm to use.
%          Non-negative values correspond to selecting the algorithm for
%          performing an LU decomposition. This corresponds to the
%          algorithm input of the LUDecompSquare function. algorithm=-1
%          means that a Cholesky-decomposition algorithm should be used.
%          The default if omitted or an empty matrix is passed is 3 unless
%          n<=5, in which case an explicit sum of terms is used. Algorithms
%          0 and 2 tend not to be as numerically stable as the others.
%          Algorithm=-1 can only be used if the matrix is known to be
%          symmetric and positive definite.
%
%OUTPUTS: val The value of the determinant.
%
%In the proof of Theorem 3.2.1 in Section 3.2.5 of [1] shows the relation
%between the product of the diagonals of U and the determinant of A. Here,
%the values on the diagonal of L are not always 1, so one must multiply by
%those as well. Also, pivoting changes the sign of the determinant, so one
%must account for the number of swaps of rows and columns due to pivots.
%
%For algorithm=-1, the squared product of the diagonals of the Cholesky
%decomposition matrix is equal to the determinant of the original matrix.
%
%EXAMPLE 1:
%Here, we show that for a 5D matrix, the explicit formula here is faster
%than Matlab's. However, one can see that this only takes effect on the
%second line calling the function, due to however, matlab's internal
%optimization works. Thus, we time Matlab's det function twice and this
%function twice. Matlab's function is the same each time. However, the
%second time this function is times, it beats Matlab's function.
% x5=rand(5,5);
% determinant(x5);
% tic%Matlab's function
% for k=1:1e6
% det(x5);
% end
% toc
% tic%Matlab's function again
% for k=1:1e6
% det(x5);
% end
% toc
% tic%This function the first time. This time will be slower than Matlab's.
% for k=1:1e6
% determinant(x5);
% end
% toc
% tic%This function again (The fastest of them all).
% for k=1:1e6
% determinant(x5);
% end
% toc
%
%EXAMPLE 2:
%Problem 4.6.1 of Chapter 4 of [1] gives an explicit formula for the
%determinant of a Vandermonde matrix. Here, we compare the relative errors
%between the explicit result and the solution with with algorithms 1, 3, 4,
%and 5 and also Matlab's built-in function when computing the determinant
%of a Vandermone matrix.
% numMCRuns=1e4;
% n=10;
% MatErr=0;
% AlgErr1=0;
% AlgErr3=0;
% AlgErr4=0;
% AlgErr5=0;
% for k=1:numMCRuns
%     alpha=randn(n,1);
%     A=vander(alpha);
%     detValTrue=1;
%     for j=1:(n-1)
%         for i=(j+1):n
%             detValTrue=detValTrue*(alpha(j)-alpha(i));
%         end
%     end
%     MatErr=MatErr+abs((detValTrue-det(A))/detValTrue);
%     detVal=determinant(A,1);
%     AlgErr1=AlgErr1+abs((detValTrue-detVal)/detValTrue);
%     detVal=determinant(A,3);
%     AlgErr3=AlgErr3+abs((detValTrue-detVal)/detValTrue);
%     detVal=determinant(A,4);
%     AlgErr4=AlgErr4+abs((detValTrue-detVal)/detValTrue);
%     detVal=determinant(A,5);
%     AlgErr5=AlgErr5+abs((detValTrue-detVal)/detValTrue);
% end
% MatErr=MatErr/numMCRuns;
% AlgErr1=AlgErr1/numMCRuns
% AlgErr3=AlgErr3/numMCRuns
% AlgErr4=AlgErr4/numMCRuns
% AlgErr5=AlgErr5/numMCRuns
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    switch(size(A,1))
        case 2
            val=A(1,1)*A(2,2)-A(1,2)*A(2,1);
            return;
        case 3
            val=-A(1,3)*A(2,2)*A(3,1)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)+A(1,1)*A(2,2)*A(3,3);
            return;
        case 4
            val=A(1,4)*A(2,3)*A(3,2)*A(4,1)-A(1,3)*A(2,4)*A(3,2)*A(4,1)-A(1,4)*A(2,2)*A(3,3)*A(4,1)+A(1,2)*A(2,4)*A(3,3)*A(4,1)+A(1,3)*A(2,2)*A(3,4)*A(4,1)-A(1,2)*A(2,3)*A(3,4)*A(4,1)-A(1,4)*A(2,3)*A(3,1)*A(4,2)+A(1,3)*A(2,4)*A(3,1)*A(4,2)+A(1,4)*A(2,1)*A(3,3)*A(4,2)-A(1,1)*A(2,4)*A(3,3)*A(4,2)-A(1,3)*A(2,1)*A(3,4)*A(4,2)+A(1,1)*A(2,3)*A(3,4)*A(4,2)+A(1,4)*A(2,2)*A(3,1)*A(4,3)-A(1,2)*A(2,4)*A(3,1)*A(4,3)-A(1,4)*A(2,1)*A(3,2)*A(4,3)+A(1,1)*A(2,4)*A(3,2)*A(4,3)+A(1,2)*A(2,1)*A(3,4)*A(4,3)-A(1,1)*A(2,2)*A(3,4)*A(4,3)-A(1,3)*A(2,2)*A(3,1)*A(4,4)+A(1,2)*A(2,3)*A(3,1)*A(4,4)+A(1,3)*A(2,1)*A(3,2)*A(4,4)-A(1,1)*A(2,3)*A(3,2)*A(4,4)-A(1,2)*A(2,1)*A(3,3)*A(4,4)+A(1,1)*A(2,2)*A(3,3)*A(4,4);
            return;
        case 5
            val=(A(1,5)*A(2,4)*A(3,3)*A(4,2)-A(1,4)*A(2,5)*A(3,3)*A(4,2)-A(1,5)*A(2,3)*A(3,4)*A(4,2)+A(1,3)*A(2,5)*A(3,4)*A(4,2)+A(1,4)*A(2,3)*A(3,5)*A(4,2)-A(1,3)*A(2,4)*A(3,5)*A(4,2)-A(1,5)*A(2,4)*A(3,2)*A(4,3)+A(1,4)*A(2,5)*A(3,2)*A(4,3)+A(1,5)*A(2,2)*A(3,4)*A(4,3)-A(1,2)*A(2,5)*A(3,4)*A(4,3)-A(1,4)*A(2,2)*A(3,5)*A(4,3)+A(1,2)*A(2,4)*A(3,5)*A(4,3)+A(1,5)*A(2,3)*A(3,2)*A(4,4)-A(1,3)*A(2,5)*A(3,2)*A(4,4)-A(1,5)*A(2,2)*A(3,3)*A(4,4)+A(1,2)*A(2,5)*A(3,3)*A(4,4)+A(1,3)*A(2,2)*A(3,5)*A(4,4)-A(1,2)*A(2,3)*A(3,5)*A(4,4)-A(1,4)*A(2,3)*A(3,2)*A(4,5)+A(1,3)*A(2,4)*A(3,2)*A(4,5)+A(1,4)*A(2,2)*A(3,3)*A(4,5)-A(1,2)*A(2,4)*A(3,3)*A(4,5)-A(1,3)*A(2,2)*A(3,4)*A(4,5)+A(1,2)*A(2,3)*A(3,4)*A(4,5))*A(5,1)-(A(1,5)*A(2,4)*A(3,3)*A(4,1)-A(1,4)*A(2,5)*A(3,3)*A(4,1)-A(1,5)*A(2,3)*A(3,4)*A(4,1)+A(1,3)*A(2,5)*A(3,4)*A(4,1)+A(1,4)*A(2,3)*A(3,5)*A(4,1)-A(1,3)*A(2,4)*A(3,5)*A(4,1)-A(1,5)*A(2,4)*A(3,1)*A(4,3)+A(1,4)*A(2,5)*A(3,1)*A(4,3)+A(1,5)*A(2,1)*A(3,4)*A(4,3)-A(1,1)*A(2,5)*A(3,4)*A(4,3)-A(1,4)*A(2,1)*A(3,5)*A(4,3)+A(1,1)*A(2,4)*A(3,5)*A(4,3)+A(1,5)*A(2,3)*A(3,1)*A(4,4)-A(1,3)*A(2,5)*A(3,1)*A(4,4)-A(1,5)*A(2,1)*A(3,3)*A(4,4)+A(1,1)*A(2,5)*A(3,3)*A(4,4)+A(1,3)*A(2,1)*A(3,5)*A(4,4)-A(1,1)*A(2,3)*A(3,5)*A(4,4)-A(1,4)*A(2,3)*A(3,1)*A(4,5)+A(1,3)*A(2,4)*A(3,1)*A(4,5)+A(1,4)*A(2,1)*A(3,3)*A(4,5)-A(1,1)*A(2,4)*A(3,3)*A(4,5)-A(1,3)*A(2,1)*A(3,4)*A(4,5)+A(1,1)*A(2,3)*A(3,4)*A(4,5))*A(5,2)+(A(1,5)*A(2,4)*A(3,2)*A(4,1)-A(1,4)*A(2,5)*A(3,2)*A(4,1)-A(1,5)*A(2,2)*A(3,4)*A(4,1)+A(1,2)*A(2,5)*A(3,4)*A(4,1)+A(1,4)*A(2,2)*A(3,5)*A(4,1)-A(1,2)*A(2,4)*A(3,5)*A(4,1)-A(1,5)*A(2,4)*A(3,1)*A(4,2)+A(1,4)*A(2,5)*A(3,1)*A(4,2)+A(1,5)*A(2,1)*A(3,4)*A(4,2)-A(1,1)*A(2,5)*A(3,4)*A(4,2)-A(1,4)*A(2,1)*A(3,5)*A(4,2)+A(1,1)*A(2,4)*A(3,5)*A(4,2)+A(1,5)*A(2,2)*A(3,1)*A(4,4)-A(1,2)*A(2,5)*A(3,1)*A(4,4)-A(1,5)*A(2,1)*A(3,2)*A(4,4)+A(1,1)*A(2,5)*A(3,2)*A(4,4)+A(1,2)*A(2,1)*A(3,5)*A(4,4)-A(1,1)*A(2,2)*A(3,5)*A(4,4)-A(1,4)*A(2,2)*A(3,1)*A(4,5)+A(1,2)*A(2,4)*A(3,1)*A(4,5)+A(1,4)*A(2,1)*A(3,2)*A(4,5)-A(1,1)*A(2,4)*A(3,2)*A(4,5)-A(1,2)*A(2,1)*A(3,4)*A(4,5)+A(1,1)*A(2,2)*A(3,4)*A(4,5))*A(5,3)-(A(1,5)*A(2,3)*A(3,2)*A(4,1)-A(1,3)*A(2,5)*A(3,2)*A(4,1)-A(1,5)*A(2,2)*A(3,3)*A(4,1)+A(1,2)*A(2,5)*A(3,3)*A(4,1)+A(1,3)*A(2,2)*A(3,5)*A(4,1)-A(1,2)*A(2,3)*A(3,5)*A(4,1)-A(1,5)*A(2,3)*A(3,1)*A(4,2)+A(1,3)*A(2,5)*A(3,1)*A(4,2)+A(1,5)*A(2,1)*A(3,3)*A(4,2)-A(1,1)*A(2,5)*A(3,3)*A(4,2)-A(1,3)*A(2,1)*A(3,5)*A(4,2)+A(1,1)*A(2,3)*A(3,5)*A(4,2)+A(1,5)*A(2,2)*A(3,1)*A(4,3)-A(1,2)*A(2,5)*A(3,1)*A(4,3)-A(1,5)*A(2,1)*A(3,2)*A(4,3)+A(1,1)*A(2,5)*A(3,2)*A(4,3)+A(1,2)*A(2,1)*A(3,5)*A(4,3)-A(1,1)*A(2,2)*A(3,5)*A(4,3)-A(1,3)*A(2,2)*A(3,1)*A(4,5)+A(1,2)*A(2,3)*A(3,1)*A(4,5)+A(1,3)*A(2,1)*A(3,2)*A(4,5)-A(1,1)*A(2,3)*A(3,2)*A(4,5)-A(1,2)*A(2,1)*A(3,3)*A(4,5)+A(1,1)*A(2,2)*A(3,3)*A(4,5))*A(5,4)+(A(1,4)*A(2,3)*A(3,2)*A(4,1)-A(1,3)*A(2,4)*A(3,2)*A(4,1)-A(1,4)*A(2,2)*A(3,3)*A(4,1)+A(1,2)*A(2,4)*A(3,3)*A(4,1)+A(1,3)*A(2,2)*A(3,4)*A(4,1)-A(1,2)*A(2,3)*A(3,4)*A(4,1)-A(1,4)*A(2,3)*A(3,1)*A(4,2)+A(1,3)*A(2,4)*A(3,1)*A(4,2)+A(1,4)*A(2,1)*A(3,3)*A(4,2)-A(1,1)*A(2,4)*A(3,3)*A(4,2)-A(1,3)*A(2,1)*A(3,4)*A(4,2)+A(1,1)*A(2,3)*A(3,4)*A(4,2)+A(1,4)*A(2,2)*A(3,1)*A(4,3)-A(1,2)*A(2,4)*A(3,1)*A(4,3)-A(1,4)*A(2,1)*A(3,2)*A(4,3)+A(1,1)*A(2,4)*A(3,2)*A(4,3)+A(1,2)*A(2,1)*A(3,4)*A(4,3)-A(1,1)*A(2,2)*A(3,4)*A(4,3)-A(1,3)*A(2,2)*A(3,1)*A(4,4)+A(1,2)*A(2,3)*A(3,1)*A(4,4)+A(1,3)*A(2,1)*A(3,2)*A(4,4)-A(1,1)*A(2,3)*A(3,2)*A(4,4)-A(1,2)*A(2,1)*A(3,3)*A(4,4)+A(1,1)*A(2,2)*A(3,3)*A(4,4))*A(5,5);
            return
        case 1
            val=A;
            return
        case 0
            val=1;
            return 
        otherwise
            algorithm=3;
    end
end

if(algorithm==-1)
    val=prod(diag(chol(A,'lower'))).^2;
else
    [L,U,~,~,sgn1,sgn2]=LUDecompSquare(A,algorithm);
    val=sgn1*sgn2*prod(diag(L))*prod(diag(U));
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
 