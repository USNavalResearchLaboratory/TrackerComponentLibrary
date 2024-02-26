function vals=chirpZTrans(x,M,W,A,algorithm)
%%CHIRPZTRANS Compute the chirp-z transform of the length-n sequence x at M
%             points. For k=0,...,M-1, this computes
%             X(k+1)=sum_{n=0}^{N-1}x(n+1)*A^(-n)*W^(n*k)
%             efficiently using FFTs. Sometimes, this sum is written as
%             X(k+1)=sum_{n=0}^{N-1}x(n+1)*zk^n
%             where
%             zk=W^k/A
%             Multiple chirp-z transforms can be computed at once by
%             passing a matrix for x. The special case of A=1, M=N and
%             W=exp(-1j*2*pi/N) is just an FFT. The chirp-z transform plays
%             a role in performing range-Doppler matched filtering
%             interpolation and in handling range migration.
%
%INPUTS: x An NXnumSeq matrix; this can be complex. The chirp-z transform
%          is performed over the columns of the matrix. If x is of type
%          gpuArray, then the transform is performed on the GPU and vals
%          will be of type GPUArray.
%        M The length of the output sequence. If an empty matrix is passed,
%          then N is used.
%        W The base for one of the exponential terms of the chirp-z
%          transform. This can be complex. This is numSeqX1 or 1XnumSeq if
%          the value is different for each column or 1X1 if it is the same
%          for all columns of x. If x is a gpuArray then this should be one
%          too.
%        A The base for another exponential term of the transform. If this
%          parameter is omitted or an empty matrix is passed, the default
%          of A=1 is used. This can be complex. This is numSeqX1 or
%          1XnumSeq if the value is different for each column or 1X1 if it
%          is the same for all columns of x. If x is a gpuArray then this
%          should be one too.
% algorithm This specified the algorithm used. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            fast algorithm of [1], which utilizes ffts.
%          1 Directly evaluate the sums from the definition of the chirp z
%            transform. This is slow.
%
%OUTPUTS: vals The MXnumSeq set of chirp-z transformed columns of x. If x
%              is a gpuArray, then this is a gpuArray.
%
%The implementation of the chirp-z transform in algorithm 0 is based on the
%algorithm of [1]. Note that the use of the GPU requires the Parallel
%Processing Toolbox in Matlab and is only supported with algorithm 0.
%
%In most applications, A and W probably have unit magnitude and algorithm 0
%will work very well. However, if A and W do not have unit magnitude, then
%it is worth noting that algorithm 0 tends to lose accuracy due to finite
%precision limits faster than algorithm 1.
%
%EXAMPLE:
%Here, we show that algorithms 0 and 1 produce the same result (judging by
%the relative error being in line with finite precision limitations) when A
%and W are unit magnitude. We then also consider the case where W and A are
%mangitude 2. In the example, algorithm 0 causes the results to blow up, so
%the values retruned by the algorithms disagree massively.
% %First, the case where  A and W have unit mangitude.
% N=20;
% M=20;
% x=randn(N,1)+1j*randn(N,1);
% A=exp(1j*2*pi*randn(1));
% W=exp(1j*2*pi*randn(1));
% X0=chirpZTrans(x,M,W,A,0);
% X1=chirpZTrans(x,M,W,A,1);
% RelErr=max(abs(X1-X0)./abs(X1))
% %Now, the same problem where W and A have magnitude 2.
% A=A*2;
% W=W*2;
% X0=chirpZTrans(x,M,W,A,0);
% X1=chirpZTrans(x,M,W,A,1);
% RelErr=max(abs(X1-X0)./abs(X1))
%
%REFERENCES:
%[1] L. R. Rabiner, R. W. Shafer, and C. M. Rader, "The chirp z transform
%    algorithm," IEEE Transactions on Audio and Electroacoustics, vol. AU-
%    17, no. 2, pp. 86-92, 1969.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=0;
end

N=size(x,1);
%The number of columns is the number of things that this is done in
%parallel.
numSeq=size(x,2);

if(isempty(M))
    M=N;
end

switch(algorithm)
    case 0
        if(nargin<3||isempty(A))
           if(isa(x,'gpuArray'))
               A=gpuArray(1);
           else
               A=1;
           end
        end

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
    case 1
        vals=zeros(M,numSeq);
        for k=0:(M-1)
            zMult=W^k/A;
            zCur=1;
            for n=0:(N-1)
                vals(k+1,:)=vals(k+1,:)+x(n+1,:)*zCur;
                zCur=zCur*zMult;
            end
        end
    otherwise
        error('Unknown algorithm selected.')
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
