function [R,t,c]=findTransParam(A,B)
%%FINDTRANSPARAM Given two sets of vectors, this algorithm find the
%                transformation parameters to transform the vectors in A to
%                B. If one output parameter is desired, the transformation
%                is just a least-squares estimate of the rotation matrix
%                for the transformation A=R*B. If more outputs are desired,
%                then a least squares estimate of the parameters for the
%                transformation A=bsxfun(@plus,c*R*B,t) is given. If only R
%                and t are desired, it is assumed that c=1.
%
%INPUTS: A An mXn matrix of n m-dimensional real vectors in the destination
%          coordinate system.
%        B An mXn matrix of n m-dimensional real vectors in the source
%          coordinate system.
%
%OUTPUTS: R A rotation matrix for the transformation from B to A.
%         t A translation vector for the transformation from B to A.
%         c A scale factor for the transformation from B to A.
%
%For both scenarios, when just R is desired and when more parameters are
%desired, the algorithms of [1] are used. When only R is desired, the
%solution in Equation 4 is used. When more parameters are desired, then the
%solution of Equation 40-42 are used. If the scale factor c is not
%requested, then the translation t is computed assuming that c is one.
%
%The cost function used in the paper minimizes the cost function 
%||A-R*B||^2, where the vertical lines indicate the sum of the L2 norm
%across columns. When a scaling and translation parameter are desired, the
%function minimizes ||A-bsxfun(@plus,c*R*B,t)||^2.
%
%This function assumes that all of the input vectors have the same
%precision. If only a rotation matrix is desired and the vectors all have
%the same magnitude (e.g. are all unit vectors), but have varying
%accuracies, the variable accuracies can be taken into account by
%pre-scaling the data. For example, if point i corresponding to vectors 
%A(:,i) and B(:,i) has an inverse covariance matrix inv(P)=S*S', then
%instead of vectors A(:,i) and B(:,i), one should use vectors A(:,i)S and
%B(:,i)S. This effectively changes the cost function from a sum of
%penalties of the form ||A(:,i)-R*B(:,i)||^2 to a sum of Mahalanobis
%distances of the form (A(:,i)-R*B(:,i))'*inv(P)*(A(:,i)-R*B(:,i)). Such a
%simple rescaling of the inputs does not work when a translation vector
%might be present.
%
%REFERENCES:
%[1] S. Umeyama, "Least squares estimation of transformation parameters
%    between two point patterns," IEEE Transactions on Pattern Analysis and
%    Machine Intelligence, vol. 13, no. 4, pp. 376-380, Apr. 1991.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    m=size(A,1);
    n=size(A,2);
    
    %If only a rotation matrix is desired.
    if(nargout<=1)
        ABp=A*B';
        [U,~,V]=svd(ABp);
        rankVal=rank(ABp);

        S=eye(m);
        if(rankVal<m-1)
            error('The matrices provided do not provide enough degrees of freedom to uniquely solve the problem.');
        elseif(rankVal==m-1)
            %Implementing Equation 5.
            if(det(U)*det(V)<0)
               S(m,m)=-1; 
            end
        else
            %Implementing Equation 3
            if(det(ABp)<0)
                S(m,m)=-1;
            end
        end

        %Equation 4.
        R=U*S*V';
    else
        muX=(1/n)*sum(B,2);%Equation 34
        muY=(1/n)*sum(A,2);%Equation 35
        diffX=bsxfun(@minus,B,muX);
        sigma2X=(1/n)*sum(sum(diffX.*diffX));%Equation 36
        diffY=bsxfun(@minus,A,muY);
        %sigma2Y=(1/n)*sum(sum(diffY.*diffY));%Equation 37
        SigmaXY=(1/n)*diffY*diffX';%Equation 38
        
        [U,D,V]=svd(SigmaXY);
        S=eye(m);
        rankVal=rank(SigmaXY);
        
        if(rankVal<m-1)
            error('The matrices provided do not provide enough degrees of freedom to uniquely solve the problem.');
        elseif(rankVal==m-1)
            %Implementing Equation 43
            if(det(U)*det(V)<0)
               S(m,m)=-1; 
            end
        else
            %Implementing Equation 39
            if(det(SigmaXY)<0)
                S(m,m)=-1;
            end
        end
        
        R=U*S*V';%Equation 40
        if(nargout==2)
            c=1;
        else
            c=(1/sigma2X)*trace(D*S);%Equation 42
        end
        t=muY-c*R*muX;%Equation 41
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
