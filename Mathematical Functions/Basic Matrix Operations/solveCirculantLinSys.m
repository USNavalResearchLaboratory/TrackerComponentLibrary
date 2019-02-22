function x=solveCirculantLinSys(z,y)
%%SOLVECIRCULANTLINSYS Solve the equation C(z)*x=y for x, where C(z) is a
%           circulant matrix where the first column is z. z and y can be
%           complex. A circulant matrix has the form
%           C(z)=[z(1),  z(n),  z(n-1), ... z(2);
%                 z(2),  z(1),  z(n),   ... z(3);
%                 ...    ...    ...     ... ...
%                 z(n-1),z(n-2),z(n-3), ... z(n);
%                 z(n),  z(n-1),z(n-2), ... z(1)];
%
%INPUTS: z The first column of the nXn circulant matrix. This can be given
%          as a row or a column vector.
%        y The nX1 right-hand side of the equation C(z)*x=y.
%
%OUTPUTS: x The solution to the system C(z)*x=y
%
%This function implements Algorithm 4.8.1 of [1], which if efficient with
%large matrices due to its use of FFTs.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

w=conj(fft(conj(y(:))))./conj(fft(conj(z(:))));
x=fft(w)/length(z);

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
