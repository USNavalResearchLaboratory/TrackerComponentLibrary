function vals=chirpZTrans(x,M,W,A)
%%CHIRPZTRANS Compute the chirp-z transform of the length-n sequence x at M
%             points. For k=0,...,M-1, this computes
%             X(k+1)=sum_{n=0}^{N-1}x(n+1)*A^(-n)*W^(n*k)
%             efficiently using FFTs. Multiple chirp-=z transforms can be
%             computed at once by passing a matrix for x. The special case
%             of A=1, M=n and W=exp(-1j*2*pi/N) is just an FFT. The chirp-z
%             transform plays a role in performing range-doppler matched
%             filtering taking into account range migration.
%
%INPUTS: x An NXnumSeq matrix; this can be complex. The chirp-z transform
%          is performed over the columns of the matrix. If x is of type
%          gpuArray, then the transform is performed on the GPU and vals
%          will be of type GPUArray.
%        M The length of the output sequence. Often, this is just N.
%        W The base for one of the exponential terms of the chirp-z
%          transform. This can be complex. This is numSeqX1 or 1XnumSeq if
%          the value is different for each column or 1X1 if it is the same
%          for all columns of x. If x is a gpuArray then this should be one
%          too.
%        A The base for another exponential term of the transform. if this
%          parameter is omitted or an empty matrix is passed, the default
%          of A=1 is used. This can be complex. This is numSeqX1 or
%          1XnumSeq if the value is different for each column or 1X1 if it
%          is the same for all columns of x. If x is a gpuArray then this
%          should be one too.
%
%OUTPUTS: vals The MXnumSeq set of chirp-z transformed columns of x. If x
%              is a gpuArray, then this is a gpuArray.
%
%The implementation of the chirp-z transform is based on the algorithm of
%[1]. Note that the use of the GPU requires the Parallel Processing Toolbox
%in Matlab. 
%
%REFERENCES:
%[1] L. R. Rabiner, R. W. Shafer, and C. M. Rader, "The chirp z transform
%    algorithm," IEEE Transactions on Audio and Electroacoustics, vol. AU-
%    17, no. 2, pp. 86-92, 1969.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(A))
       if(isa(x,'gpuArray'))
           A=gpuArray(1);
       else
           A=1;
       end
    end

    %The number of transformed points
    N=size(x,1);
    %The number of columns is the number of things that this is done on in
    %parallel.
    numSeq=size(x,2);

    %Round L to the next highest power of two.
    L=N+M-1;
    [F,E]=log2(L);
    %If it is not already a power of two, round up to the next power of
    %two.
    if(F~=0.5)
        L=2^E;
    end

    if(isscalar(W))
       W=repmat(W,1,numSeq);
    else
       W=reshape(W,1,numSeq);
    end
    
    if(isscalar(A))
       A=repmat(A,1,numSeq);
    else
       A=reshape(A,1,numSeq);
    end

    if(isa(x,'gpuArray')==false)
        %Do everything on the CPU
        y=zeros(L,numSeq);
        n=(0:(N-1)).';
        y(1:N,:)=bsxfun(@power,A,-n).*bsxfun(@power,W,(n.^2/2)).*x;
        clear n x
        y=fft(y,[],1);

        %The arbitrary part of v is just set to zero.
        v=zeros(L,numSeq);
        n=(0:(M-1)).';
        v(n+1,:)=bsxfun(@power,W,-n.^2/2);
        n=((L-N+1):(L-1)).';
        v(n+1,:)=bsxfun(@power,W,-(L-n).^2/2);
        clear n
        v=fft(v,[],1);

        g=y.*v;
        clear y v
        g=ifft(g,[],1);

        g=g(1:M,:);%The extra terms are discarded.
        n=((0:(M-1))).';
        vals=bsxfun(@power,W,n.^2/2).*g;
        clear W g n
    else%Do the transformation on the GPU
        y=gpuArray.zeros(L,numSeq);
        n=gpuArray.colon(0,(N-1)).';
        y(1:N,:)=bsxfun(@power,A,-n).*bsxfun(@power,W,(n.^2/2)).*x;
        clear n x
        y=fft(y,[],1);

        %The arbitrary part of v is just set to zero.
        v=gpuArray.zeros(L,numSeq);
        n=gpuArray.colon(0,(M-1)).';
        v(n+1,:)=bsxfun(@power,W,-n.^2/2);
        n=gpuArray.colon((L-N+1),(L-1)).';
        v(n+1,:)=bsxfun(@power,W,-(L-n).^2/2);
        clear n
        v=fft(v,[],1);

        g=y.*v;
        clear y v
        g=ifft(g,[],1);

        g=g(1:M,:);%The extra terms are discarded.
        n=gpuArray.colon(0,(M-1)).';
        vals=bsxfun(@power,W,n.^2/2).*g;
        clear W g n
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
