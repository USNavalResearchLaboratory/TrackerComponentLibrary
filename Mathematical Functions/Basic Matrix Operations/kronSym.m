function K=kronSym(A,B)
%%KRONSYM Take the symmetric Kronecker product of the matrices A and B. The
%         standard Kronecker product of two real, square matrices A and B,
%         one can write the following relation with respect to any real
%         square matrix S:
%         kron(A,B)*vec(S)==vec(B*S*A')
%         However, the symmetric Kronecker product is such that for any two
%         real square matrices A and B and a square SYMMETRIC matrix S, one
%         can write
%         kronSym(A,B)*vech(S,sqrt(2))==vech((1/2)*(B*S*A'+A*S*B'),sqrt(2))
%         Whereas is A and B are nXn, the standard kronecker product is
%         n^2Xn^2, the symmetric Kronecker product is
%         (n*(n+1)/2)X(n*(n+1)/2).
%
%INPUTS: A, B Two real nXn matrices. They need not be square.
%
%OUTPUTS: K The (n*(n+1)/2)X(n*(n+1)/2) symmetric Kronecker product of A
%           and B.
%
%Symmetric Kronecker products are discussed in the appendix of [1], where
%they play a role in the implementation of a semidefinite programming
%algorithm.
%
%The symmetric Kronecker product is implemented by noting that in the
%identity
%kronSym(A,B)*vech(S,sqrt(2))==vech((1/2)*(B*S*A'+A*S*B'),sqrt(2)), if S is
%a symmetric matrix where on or below the main diagonal there is only a
%single nonzero element, then one effectively selects a column of
%kronSym(A,B). To get an unscaled column of kronSym(A,B), the nonzero
%element should be 1 if it is on the main diagonal and it should be
%1/sqrt(2) for elements below (and above) the main diagonal due to the
%sqrt(2) scaling in the vech command. Basically, a simple implementation
%would be 
% n=size(A,1);
% prodDim=n*(n+1)/2;
% 
% K=zeros(prodDim,prodDim);
% SelMat=zeros(n,n);
% for curEl=1:prodDim
%     [row,col]=vechInd2Sub(n,curEl);
%     if(row==col)
%         SelMat(row,col)=1;
%     else%row>col
%         SelMat(row,col)=1/sqrt(2);
%         SelMat(col,row)=1/sqrt(2);
%     end
%     K(:,curEl)=vech((1/2)*(B*SelMat*A'+A*SelMat*B'),sqrt(2));
% 
%     SelMat(row,col)=0;
%     SelMat(col,row)=0;
% end
%However, that is rather slow. Thus, this implementation tries to directly
%write the elements of the output matrix K. 
%
%REFERENCES:
%[1] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
%    interior-point methods for semidefinite programming: Convergence
%    rates, stability and numerical results," SIAM Journal on Optimization,
%    vol. 8, no. 3, pp. 746-768, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The columns of K are indexed using (i1,i2) and the rows are indexed using
%(j1,j2).
n=size(A,1);
prodDim=n*(n+1)/2;

K=zeros(prodDim,prodDim);
sqrt2=sqrt(2);
dblSqrt2=2*sqrt2;
for i1=1:n
    col=i1+(i1-1)*(2*n-i1)/2;
    for j1=1:n
        row=j1+(j1-1)*(2*n-j1)/2;
        
        K(row,col)=(A(j1,i1)*B(j1,i1)+B(j1,i1)*A(j1,i1))/2;
        for j2=(j1+1):n
            row=row+1;
            
            K(row,col)=(A(j1,i1)*B(j2,i1)+B(j1,i1)*A(j2,i1))/sqrt2;
        end
    end
    
    for i2=(i1+1):n
      col=col+1;
      for j1=1:n
         row=j1+(j1-1)*(2*n-j1)/2;
         
         K(row,col)=(A(j1,i1)*B(j1,i2)+A(j1,i2)*B(j1,i1)+B(j1,i1)*A(j1,i2)+B(j1,i2)*A(j1,i1))/dblSqrt2;     
         for j2=(j1+1):n
            row=row+1;
            
            K(row,col)=(A(j1,i1)*B(j2,i2)+A(j1,i2)*B(j2,i1)+B(j1,i1)*A(j2,i2)+B(j1,i2)*A(j2,i1))/2;
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
