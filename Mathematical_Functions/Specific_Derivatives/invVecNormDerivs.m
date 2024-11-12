function [J,H,M]=invVecNormDerivs(cVec,dc,d2c,d3c)
%%INVVECNORMDERIVS Consider the inverse of the norm of a vector
% 1/norm(cVec). The value c can be a function of one or more variables.
% This function evaluates the first (and up to the third if requested)
% derivatives of the function with respect to the values on which cVec
% depends.
%
%INPUTS: cVec An NX1 vector.
%          dc An NXm set of NX1 vectors where dc(:,i) holds the partial
%             derivative of cVec with respect to the ith parameter.
%         d2c An NXmXm set of NX1 vectors where d2c(:,i,j) holds the
%             second derivative of cVec with respect to the ith and jth
%             parameters (which can be the same). This is only used if more
%             than one output is requested.
%         d3c An NXmXmXm set of NX1 vectors where d3c(:,i,j,k) holds the
%             third derivative of cVec with respect to the ith, jth and kth
%             parameters (which can be unique or not). This is only used if
%             more than two outputs are requested.
%           
%OUTPUTS: J A 1Xm set of first derivatives of 1/norm(cVec), with J(i)
%           holding the derivative with respect to the ith parameter.
%         H An mXm set of second derivatives of 1/norm(cVec) with H(i,j)
%           holding the second derivative with respect to parameters i and
%           j.
%         M An mXmXm set of third derivatives of 1/norm(cVec) with M(i,j,k)
%           holding the third derivatives with respect to parameters i, j,
%            and k.
%
%This computes the value of the derivatives with an explicit solution.
%
%EXAMPLE 1:
%In this example, the first derivatives returned by this function are
%compared to finite differencing for the case where a 2D vector vectors is
%a nonlinear functions of two variables. The results agree to more than 7
%digits, which is good.
% epsVal=1e-8;
% ab=[-4;5];
% a=ab(1);
% b=ab(2);
% 
% c=[a*b;
%     b^2*a^3];
% dc=[        b, a;
%      3*b^2*a^2,2*b*a^3];
% J=invVecNormDerivs(c,dc);
% 
% invNorm=1/norm(c);
% invNormEps=zeros(1,2);
% for curEpsVar=1:2
%     abCur=ab;
%     abCur(curEpsVar)=abCur(curEpsVar)+epsVal;
%     a=abCur(1);
%     b=abCur(2);
% 
%     cEps=[a*b;
%            b^2*a^3];
%     invNormEps(curEpsVar)=1/norm(cEps);
% end
% JNumDiff=zeros(1,2);
% for k=1:2
%     JNumDiff(k)=(invNormEps(k)-invNorm)/epsVal;
% end
% RelErr=max(abs((JNumDiff-J)./JNumDiff))
%
%EXAMPLE 2:
%In this example, the second derivatives returned by this function are
%compared to finite differencing for the case where a 2D vector vectors is
%a nonlinear functions of two variables. The results agree to more than 7
%digits, which is good.
% epsVal=1e-8;
% ab=[-4;5];
% a=ab(1);
% b=ab(2);
% 
% c=[a*b;
%     b^2*a^3];
% dc=[        b, a;
%      3*b^2*a^2,2*b*a^3];
% d2c=zeros(2,2,2);
% d2c(1,1,:)=[0;6*b^2*a];
% d2c(1,2,:)=[1;6*b*a^2];
% d2c(2,1,:)=d2c(1,2,:);
% d2c(2,2,:)=[0;2*a^3];
% [J,H]=invVecNormDerivs(c,dc,d2c);
% 
% JEps=zeros(1,2,2);
% for curEpsVar=1:2
%     abCur=ab;
%     abCur(curEpsVar)=abCur(curEpsVar)+epsVal;
%     a=abCur(1);
%     b=abCur(2);
% 
%     c=[a*b;
%        b^2*a^3];
%     dc=[        b, a;
%         3*b^2*a^2,2*b*a^3];
%     JEps(1,:,curEpsVar)=invVecNormDerivs(c,dc);
% end
% HNumDiff=zeros(2,2);
% for k1=1:2
%     JCur=J(k1);
%     for k2=1:2
%         JEpsCur=JEps(1,k1,k2);
%         HNumDiff(k1,k2)=(JEpsCur-JCur)/epsVal;
%     end
% end
% RelErr=max(abs((HNumDiff(:)-H(:))./HNumDiff(:)))
%
%EXAMPLE 3:
%In this example, the third derivatives returned by this function are
%compared to finite differencing for the case where a 2D vector vectors is
%a nonlinear functions of two variables. The results agree to more than 7
%digits, which is good.
% epsVal=1e-8;
% ab=[-4;5];
% a=ab(1);
% b=ab(2);
% 
% c=[a*b;
%     b^2*a^3];
% dc=[        b, a;
%      3*b^2*a^2,2*b*a^3];
% d2c=zeros(2,2,2);
% d2c(1,1,:)=[0;6*b^2*a];
% d2c(1,2,:)=[1;6*b*a^2];
% d2c(2,1,:)=d2c(1,2,:);
% d2c(2,2,:)=[0;2*a^3];
% %Third derivatives
% d3c=zeros(2,2,2,2);
% d3c(1,1,1,:)=[0;6*b^2];
% d3c(1,1,2,:)=[0;12*b*a];
% d3c(1,2,2,:)=[0;6*a^2];
% d3c(2,2,2,:)=[0;0];
% %All of the permutations of the indices are the same.
% d3c(2,1,1,:)=d3c(1,1,2,:);
% d3c(1,2,1,:)=d3c(1,1,2,:);
% d3c(2,1,2,:)=d3c(1,2,2,:);
% d3c(2,2,1,:)=d3c(1,2,2,:);
% [~,H,M]=invVecNormDerivs(c,dc,d2c,d3c);
% 
% HEps=zeros(2,2,2);
% for curEpsVar=1:2
%     abCur=ab;
%     abCur(curEpsVar)=abCur(curEpsVar)+epsVal;
%     a=abCur(1);
%     b=abCur(2);
% 
%     c=[a*b;
%        b^2*a^3];
%     dc=[        b, a;
%         3*b^2*a^2,2*b*a^3];
%     d2c(1,1,:)=[0;6*b^2*a];
%     d2c(1,2,:)=[1;6*b*a^2];
%     d2c(2,1,:)=d2c(1,2,:);
%     d2c(2,2,:)=[0;2*a^3];
%     [~,HEps(:,:,curEpsVar)]=invVecNormDerivs(c,dc,d2c);
% end
% MNumDiff=zeros(2,2,2);
% for k1=1:2
%     for k2=1:2
%         HCur=H(k1,k2);
%         for k3=1:2
%             HEpsCur=HEps(k1,k2,k3);
%             MNumDiff(k1,k2,k3)=(HEpsCur-HCur)/epsVal;
%         end
%     end
% end
% RelErr=max(abs((MNumDiff(:)-M(:))./MNumDiff(:)))
%
%July 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDerivVars=size(dc,2);

c2=cVec'*cVec;

J=zeros(1,numDerivVars);
for a=1:numDerivVars
    J(a)=-(cVec'*dc(:,a))/c2^(3/2);
end

if(nargout>1)
    H=zeros(numDerivVars,numDerivVars);
    for a=1:numDerivVars
        ca=dc(:,a);
    
        for b=a:numDerivVars
            cb=dc(:,b);
            cab=vec(d2c(a,b,:));
    
            H(a,b)=(3/((c2)^(5/2))*((cVec'*cb)*(cVec'*ca))-1/(c2)^(3/2)*((ca'*cb)+(cVec'*cab)));

            H(b,a)=H(a,b);
        end
    end

    if(nargout>2)
        M=zeros(numDerivVars,numDerivVars,numDerivVars);
        for a=1:numDerivVars
            ca=dc(:,a);
        
            for b=a:numDerivVars
                cb=dc(:,b);
                cab=vec(d2c(a,b,:));
        
                for c=b:numDerivVars
                    cc=dc(:,c);
                    cac=vec(d2c(a,c,:));
                    cbc=vec(d2c(b,c,:));
                    cabc=vec(d3c(a,b,c,:));

                    M(a,b,c)=-15/(c2)^(7/2)*((cVec'*cc)*(cVec'*cb)*(cVec'*ca))...
                         +3/(c2)^(5/2)*(((cb'*cc)+(cVec'*cbc))*(cVec'*ca)+(cVec'*cb)*((ca'*cc)+(cVec'*cac))+(cVec'*cc)*((ca'*cb)+(cVec'*cab)))...
                         -1/(c2)^(3/2)*(((ca'*cbc)+(cb'*cac)+(cc'*cab)+(cVec'*cabc)));
                    %All of the permutations of the indices are the same.
                    M(b,a,c)=M(a,b,c);
                    M(c,a,b)=M(a,b,c);
                    M(a,c,b)=M(a,b,c);
                    M(b,c,a)=M(a,b,c);
                    M(c,b,a)=M(a,b,c);
                end
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
