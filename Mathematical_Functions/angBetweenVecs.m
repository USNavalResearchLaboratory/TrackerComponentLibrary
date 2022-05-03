function angDiff=angBetweenVecs(v1,v2)
%%ANGBETWEENVECS Given vector pairs in 2D or higher space, find the angular
%                distances between them (non-negative values). 
%
%INPUTS: v1 A 3XN matrix of N real vectors or a single 3X1 vector is they
%           are all the same and v2 varies.
%        v2 A 3XN matrix of N real vectors or a single 3X1 vector if all of
%           them are the same and v1 varies.
%
%OUTPUTS: angDiff An 1XN matrix of the angular differences (in radians)
%                 between the vectors in v1 and the correspinding vectors
%                 in v2. The distances can range from 0 to pi.
%
%In 3D, a cross-product formula is used instead of the more logical dot-
%product formula so as to maximize numerical precision when vectors are
%nearly parallel. The vectors do not need to have the same magnitudes.
%Similarly, in 2D, a tangent formula is used. In more than 3D, a cross
%product formula is used.
%
%If one or both of the vectors has zero magnitude, then the angular
%difference will be zero. Note that angBetweenVecs(v1,v2) is the same as
%angBetweenVecs(v2,v1).
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of vectors.
N1=size(v1,2);
N2=size(v2,2);
N=max(N1,N2);
if(N1~=N)
    v1=repmat(v1,[1,N]); 
elseif(N2~=N)
    v2=repmat(v2,[1,N]); 
end

numDim=size(v1,1);

angDiff=zeros(1,N);
if(numDim==3)
    for curVec=1:N
        angDiff(curVec)=atan2(norm(cross(v1(:,curVec),v2(:,curVec))),dot(v1(:,curVec),v2(:,curVec)));
    end
elseif(numDim==2)
    for curVec=1:N
        a1=v1(1,curVec);
        a2=v1(2,curVec);
        b1=v2(1,curVec);
        b2=v2(2,curVec);

        angDiff(curVec)=atan2(abs(b2*a1-b1*a2),b1*a1+b2*a2);
    end
else
    for curVec=1:N
        norm1=norm(v1(:,curVec));
        norm2=norm(v2(:,curVec));
        
        if(norm1==0||norm2==0)
            angDiff(curVec)=0;
        else
            argVal=dot(v1(:,curVec),v2(:,curVec))/(norm1*norm2);
            %Clip to the vlaid range to deal with finite precision issues.
            argVal=max(min(argVal,1),-1);
            angDiff(curVec)=acos(argVal);
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
