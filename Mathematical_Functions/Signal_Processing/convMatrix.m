function X=convMatrix(x,N,shape)
%%CONVMATRIX Create a matrix X from the vector x such that X*y, where y is
%            a length-N vector, is the same as conv(x,y). That is, it
%            performs the linear convolution. Vectors are treated as column
%            vectors. Shape options similar to those in conv are available.
%            The linear convolution of vectors x and y is defined as
%            w(k)=sum_{j} x(j)*y(k-j+1)
%            where values for indices outside of the length of x or y are
%            taken to be zero. The sum is over all nonzero values for a
%            given k.
%
%INPUTS: x A 1XM or MX1 vector. This is treated as a column vector when
%          computing the matrix.
%        N The length of the vector with which x is convolved.
%    shape An optional parameter that determines the shape of the output
%          matrix x. Possible values are
%          'full' (The default if omitted or an empty matrix is passed).
%                 Return the full convolution. X has N+M-1 rows.
%          'same' Return only the length M central portion of the
%                 convolution. Thus X has M rows. The first row corresponds
%                 to row number fix(N/2)+1 of the full matrix.
%         'valid' Only return the part of the convolution that has no zero
%                 padding. This corresponds to rows N through M of the full
%                 matrix if N>M, then an empty matrix is returned for X.
%
%OUTPUTS: X The convolution matrix such that X*y, where y is a length-N
%          column vector, is equal to conv(x(:),y,shape).
%
%Linear convolutions are extensively used in basic signal processing
%applications and a discussion of them can be found in introductory
%textbooks.
%
%EXAMPLE:
%Here, we just show that the result is equivalent to the conv function.
% x=(-10:10)';
% y=(7:18)';
% N=length(y);
% convMatrix(x,N)*y-conv(x,y)
%One will see that the difference is all zeros.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(shape))
    shape='full'; 
end

x=x(:);
X=toeplitz([x;zeros(N-1,1)],[x(1);zeros(N-1,1)]);

switch(shape)
    case 'full'
    case 'same'
        M=length(x);
        
        %Take only the central portion
        startIdx=fix(N/2)+1;
        X=X(startIdx:(startIdx+M-1),:);
    case 'valid'
        %Remove the zero-padded edges.
        M=length(x);
        X=X(N:M,:);
    otherwise
        error('Unknown shape option specified')     
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
