function FMats=getTransMats(F,baseIdx)
%%GETTRANSMATS Given a set of discrete-time state transition matrices,
%              obtain a hypermatrix of transition matrices to prediction
%              forward and backward in time between arbitrary points.
%
%INPUTS: F An xDimXxDimXN matrix of transition matrices across the times
%          where there are measurements.
%  baseIdx As opposed to returning all of the transition matrices, if the
%          baseIdx parameter is given, then of the complete FMats matrix,
%          only FMats(:,:,:,baseIdx) is returned. 
%
%OUTPUTS: FMats An xDimXxDimXNXN hypermatrix of the F transition values,
%               unless baseIdx is given in which case it is an xDimXxDimXN
%               matrix.
%
%If x(:,b) is the target state at discrete-time step b, then, in the
%absence of process noise, the predicted state at discrete time-step a is
%given by
%FMats(:,:,a,b)*x(:,b)
%Using the baseIdx parameter set to b, one would use instead
%FMats(:,:,a)*x(:,b)
%
%The transition matrices are derived as the Phi function in [1].
%
%REFERENCES:
%[1] A. B. Poore, B. J. Slocumb, B. J. Suchomel, F. H. Obermeyer, S. M.
%    Herman, and S. M. Gadaleta, "Batch maximum likelihood (ML) and maximum
%    a posteriori (MAP) estimation with process noise for tracking
%    applications," in Proceedings of SPIE: Signal and Data Processing of
%    Small Targets, vol. 5204, San Diego, CA, 3 Aug. 2003, pp. 188-199.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(F,1);
    N=size(F,3)+1;

    FMats=zeros(xDim,xDim,N,N);

    if(nargin<2)
        baseIdx=[];
    end
    
    if(isempty(baseIdx))
        %First, fill in the diagonal matrices
        for n=1:N
            FMats(:,:,n,n)=eye(xDim);
        end

        %Next, fill in the matrices where p>i
        for i=1:(N-1)
            FRecur=F(:,:,i);
            for p=(i+1):N
                FMats(:,:,p,i)=FRecur;
                if(p<N)
                    FRecur=F(:,:,p)*FRecur;
                end
            end
        end

        %Finally, fill in the matrices where p<i
        for p=1:(N-1)
            FRecur=F(:,:,p);
            for i=(p+1):N
                FMats(:,:,p,i)=inv(FRecur);
                if(i<N)
                    FRecur=F(:,:,i)*FRecur;
                end
            end
        end
    else
        FMats=zeros(xDim,xDim,N);

        %The case where p=baseIndex
        FMats(:,:,baseIdx)=eye(xDim);
        
        %The case where p>baseIndex
        if(baseIdx<N)
            FRecur=F(:,:,baseIdx);
        end
        for p=(baseIdx+1):N
            FMats(:,:,p)=FRecur;
            if(p<N)
                FRecur=F(:,:,p)*FRecur;
            end
        end

        %The case where p<baseIndex
        for p=1:(baseIdx-1)
            FRecur=F(:,:,p);
            for i=(p+1):baseIdx
                if(i==baseIdx)
                    FMats(:,:,p)=inv(FRecur);
                end
                if(i<N)
                    FRecur=F(:,:,i)*FRecur;
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
