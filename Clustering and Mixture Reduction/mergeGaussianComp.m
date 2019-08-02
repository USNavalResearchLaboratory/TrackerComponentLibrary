function [w,mu,P]=mergeGaussianComp(w,mu,P,clust)
%%MERGEGAUSSIANCOMP Given a Gaussian mixture, replace the selected
%                   components with their weighted mean and covariance
%                   matrix.
%
%INPUTS: w A NX1 vector of weights of the components of the Gaussian
%          mixture.
%       mu An xDimXN matrix of the means of the vector components of the
%          Gaussian mixture.
%        P An xDim XxDim XN hypermatrix of the covariance matrices for the
%          components of the Gaussian mixture.
%    clust A vector of the indices of the components in the mixture that
%          are to be merged. If this parameter is omitted, then it is
%          assumed that all of the components should be merged.
%
%OUTPUTS: w The weights of the mixture after merging the selected
%           components.
%        mu The means of the mixture after merging the selected components.
%         P The covariance matrices of the mixture after merging the
%           selected components.
%
%An expression for the moments of a mixture is given in Chapter 1 of [1].
%The simplified form of the covariance sum in the book is not used, since
%precision problems can cause the merged covariance matrices to fail to
%remain positive definite. Rather, a less efficient quadratic formulation
%is used.
%
%When given a subset of components of a mixture to merge, the components
%being merged are treated as a submixture on their own. That is, they are
%merged like a normal mixture, but their weights are normalized for the
%merger. The weight of the merged component in the new mixture is the sum
%of the weights of the components being merged.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    N=length(w);
    
    if(nargin<4)
        clust=1:N;
    end
    
    numToMerge=length(w(clust));
    
    wMerged=sum(w(clust));
    wBar=w(clust)/wMerged;
    
    %Deal with numerical problems.
    if(sum(~isfinite(wBar))~=0)
       wBar(:)=1/length(wBar);
    end
    
    muMerged=calcMixtureMoments(mu(:,clust),wBar);
    
    muToMerge=mu(:,clust);
    PToMerge=P(:,:,clust);
    PMerged=0;
    for curItem=1:numToMerge
        diff=muToMerge(:,curItem)-muMerged;
        PMerged=PMerged+wBar(curItem)*(PToMerge(:,:,curItem)+diff*diff');
    end
    
    %Remove the elements that are being merged.
    sel=true(N,1);
    sel(clust)=false;
    w=w(sel);
    mu=mu(:,sel);
    P=P(:,:,sel);
    
    %Add the merged result to the end of the data set.
    w(end+1)=wMerged;
    mu(:,end+1)=muMerged;
    P(:,:,end+1)=PMerged;
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
