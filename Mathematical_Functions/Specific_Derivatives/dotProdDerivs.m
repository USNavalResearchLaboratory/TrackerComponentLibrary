function [J,H,M]=dotProdDerivs(u1,u2,du1,du2,d2u1,d2u2,d3u1,d3u2)
%%DOTPRODDERIVS Assuming that the mX1 vectors u1 and u2 are functions of
%   some other variables, this function can evaluate the first through
%   and up to the third partial derivatives of dot(u1,u2) given partial
%   derivatives of u1 and u2.
%
%INPUTS: u1 A 3X1 vector.
%        u2 A 3X1 vector.
%   du1,du2 A 3XN sets of partial first derivatives of u1 and u2,
%           respectively, with respect to N independent variables 
%           (specified by column).
%  d2u1, d2u2 An NXNX3 set of partial second derivatives of u1 and u2,
%           respectively, where d2u1(i,j,:) holds the derivative of the
%           vector with respect to the ith and jth independent variables. 
%  d3u1, d3u2 An NXNX3 set of partial third derivatives of u1 and u2,
%           respectively, where d3u1(i,j,k,:) holds the derivative of the
%           vector with respect to the ith, jth  and kth independent
%           variables. 
%
%OUTPUTS: J A 1XN set of first derivatives of dot(u1,u2), with J(i)
%           holding the derivative with respect to the ith parameter.
%         H An NXN set of second derivatives of dot(u1,u2) with H(i,j)
%           holding the second derivative with respect to parameters i and
%           j.
%         M An NXNXN set of third derivatives of dot(u1,u2) with M(i,j,k)
%           holding the third derivatives with respect to parameters i, j,
%            and k.
%
%This function just returns analytic first second and third deriavtives of
%dot(u,1,u2), where u1 and u2 depend on some variables.
%
%EXAMPLE 1:
%In this example, the first derivatives returned by this function are
%compared to finite differencing for the case where the two vectors are
%nonlinear functions of two variables. The results agree to more than 7
%digits, which is good.
% epsVal=1e-8;
% ab=[-4;5];
% a=ab(1);
% b=ab(2);
% 
% u1=[a*b;
%     b^2*a^3];
% u2=[b^3*a;
%     a^3];
% du1=[        b, a;
%      3*b^2*a^2,2*b*a^3];
% du2=[b^3,  3*b^2*a;
%    3*a^2,        0];
% J=dotProdDerivs(u1,u2,du1,du2);
% 
% uDot=dot(u1,u2);
% uDotEps=zeros(1,2);
% for curEpsVar=1:2
%     abCur=ab;
%     abCur(curEpsVar)=abCur(curEpsVar)+epsVal;
%     a=abCur(1);
%     b=abCur(2);
% 
%     u1Eps=[a*b;
%            b^2*a^3];
%     u2Eps=[b^3*a;
%            a^3];
%     uDotEps(curEpsVar)=dot(u1Eps,u2Eps);
% end
% JNumDiff=zeros(1,2);
% for k=1:2
%     JNumDiff(k)=(uDotEps(k)-uDot)/epsVal;
% end
% RelErr=max(abs((JNumDiff-J)./JNumDiff))
%
%EXAMPLE 2:
%In this example, the second derivatives returned by this function are
%compared to finite differencing for the case where the two vectors are
%nonlinear functions of two variables. The results agree to more than 7
%digits, which is good.
% epsVal=1e-8;
% ab=[-4;5];
% a=ab(1);
% b=ab(2);
% 
% u1=[a*b;
%     b^2*a^3];
% u2=[b^3*a;
%     a^3];
% du1=[        b, a;
%      3*b^2*a^2,2*b*a^3];
% du2=[b^3,  3*b^2*a;
%    3*a^2,        0];
% d2u1=zeros(2,2,2);
% d2u1(1,1,:)=[0;6*b^2*a];
% d2u1(1,2,:)=[1;6*b*a^2];
% d2u1(2,1,:)=d2u1(1,2,:);
% d2u1(2,2,:)=[0;2*a^3];
% d2u2=zeros(2,2,2);
% d2u2(1,1,:)=[0;6*a];
% d2u2(1,2,:)=[3*b^2;0];
% d2u2(2,1,:)=d2u2(1,2,:);
% d2u2(2,2,:)=[6*b*a;0];
% [J,H]=dotProdDerivs(u1,u2,du1,du2,d2u1,d2u2);
% 
% JDotEps=zeros(1,2,2);
% for curEpsVar=1:2
%     abCur=ab;
%     abCur(curEpsVar)=abCur(curEpsVar)+epsVal;
%     a=abCur(1);
%     b=abCur(2);
% 
%     u1=[a*b;
%         b^2*a^3];
%     u2=[b^3*a;
%         a^3];
%     du1=[        b, a;
%          3*b^2*a^2,2*b*a^3];
%     du2=[b^3,  3*b^2*a;
%        3*a^2,        0];
%     JDotEps(:,:,curEpsVar)=dotProdDerivs(u1,u2,du1,du2);
% end
% HNumDiff=zeros(2,2);
% for k1=1:2
%     JCur=J(:,k1);
% 
%     for k2=1:2
%         JEps=JDotEps(:,k1,k2);
%         HNumDiff(k1,k2)=(JEps-JCur)/epsVal;
%     end
% end
% RelErr=max(abs((HNumDiff(:)-H(:))./HNumDiff(:)))
%
%EXAMPLE 3:
%In this example, the third derivatives returned by this function are
%compared to finite differencing for the case where the two vectors are
%nonlinear functions of two variables. The results agree to more than 7
%digits, which is good.
% epsVal=1e-8;
% ab=[-4;5];
% a=ab(1);
% b=ab(2);
% 
% u1=[a*b;
%     b^2*a^3];
% u2=[b^3*a;
%     a^3];
% du1=[        b, a;
%      3*b^2*a^2,2*b*a^3];
% du2=[b^3,  3*b^2*a;
%    3*a^2,        0];
% d2u1=zeros(2,2,2);
% d2u1(1,1,:)=[0;6*b^2*a];
% d2u1(1,2,:)=[1;6*b*a^2];
% d2u1(2,1,:)=d2u1(1,2,:);
% d2u1(2,2,:)=[0;2*a^3];
% d2u2=zeros(2,2,2);
% d2u2(1,1,:)=[0;6*a];
% d2u2(1,2,:)=[3*b^2;0];
% d2u2(2,1,:)=d2u2(1,2,:);
% d2u2(2,2,:)=[6*b*a;0];
% 
% %Third derivatives
% d3u1=zeros(2,2,2,2);
% d3u1(1,1,1,:)=[0;6*b^2];
% d3u1(1,1,2,:)=[0;12*b*a];
% d3u1(1,2,2,:)=[0;6*a^2];
% d3u1(2,2,2,:)=[0;0];
% d3u2=zeros(2,2,2,2);
% d3u2(1,1,1,:)=[0;6];
% d3u2(1,1,2,:)=[0;0];
% d3u2(1,2,2,:)=[6*b;0];
% d3u2(2,2,2,:)=[6*a;0];
% %All of the permutations of the indices are the same.
% d3u1(1,1,2,:)=d3u1(1,1,2,:);
% d3u1(1,2,1,:)=d3u1(1,1,2,:);
% d3u1(2,1,2,:)=d3u1(1,1,2,:);
% d3u1(1,2,2,:)=d3u1(1,2,2,:);
% d3u1(2,1,2,:)=d3u1(1,2,2,:);
% d3u1(2,2,1,:)=d3u1(1,2,2,:);
% d3u2(1,1,2,:)=d3u2(1,1,2,:);
% d3u2(1,2,1,:)=d3u2(1,1,2,:);
% d3u2(2,1,2,:)=d3u2(1,1,2,:);
% d3u2(1,2,2,:)=d3u2(1,2,2,:);
% d3u2(2,1,2,:)=d3u2(1,2,2,:);
% d3u2(2,2,1,:)=d3u2(1,2,2,:);
% [~,H,M]=dotProdDerivs(u1,u2,du1,du2,d2u1,d2u2,d3u1,d3u2);
% 
% HDotEps=zeros(2,2,2);
% for curEpsVar=1:2
%     abCur=ab;
%     abCur(curEpsVar)=abCur(curEpsVar)+epsVal;
%     a=abCur(1);
%     b=abCur(2);
% 
%     u1=[a*b;
%         b^2*a^3];
%     u2=[b^3*a;
%         a^3];
%     du1=[        b, a;
%          3*b^2*a^2,2*b*a^3];
%     du2=[b^3,  3*b^2*a;
%        3*a^2,        0];
%     d2u1(1,1,:)=[0;6*b^2*a];
%     d2u1(1,2,:)=[1;6*b*a^2];
%     d2u1(2,1,:)=d2u1(1,2,:);
%     d2u1(2,2,:)=[0;2*a^3];
%     d2u2(1,1,:)=[0;6*a];
%     d2u2(1,2,:)=[3*b^2;0];
%     d2u2(2,1,:)=d2u2(1,2,:);
%     d2u2(2,2,:)=[6*b*a;0];
%     [~,HDotEps(:,:,curEpsVar)]=dotProdDerivs(u1,u2,du1,du2,d2u1,d2u2);
% end
% MNumDiff=zeros(2,2,2);
% for k1=1:2
%     for k2=1:2
%         HCur=H(k1,k2);
%         for k3=1:2
%             HEps=HDotEps(k1,k2,k3);
%             MNumDiff(k1,k2,k3)=(HEps-HCur)/epsVal;
%         end
%     end
% end
% RelErr=max(abs((MNumDiff(:)-M(:))./MNumDiff(:)))
%
%July 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDerivVars=size(du1,2);

J=zeros(1,numDerivVars);
for a=1:numDerivVars
    u1a=du1(:,a);
    u2a=du2(:,a);

    J(a)=dot(u1a,u2)+dot(u1,u2a);
end

if(nargout>1)
    H=zeros(numDerivVars,numDerivVars);
    for a=1:numDerivVars
        u1a=du1(:,a);
        u2a=du2(:,a);
    
        for b=a:numDerivVars
            u1b=du1(:,b);
            u2b=du2(:,b);
            u1ab=vec(d2u1(a,b,:));
            u2ab=vec(d2u2(a,b,:));
    
            H(a,b)=dot(u1ab,u2)+dot(u1a,u2b)+dot(u1b,u2a)+dot(u1,u2ab);
    
            %Using symmetry.
            H(b,a)=H(a,b);
        end
    end

    if(nargout>2)
        M=zeros(numDerivVars,numDerivVars,numDerivVars);
        for a=1:numDerivVars
            u1a=du1(:,a);
            u2a=du2(:,a);
    
            for b=a:numDerivVars
                u1b=du1(:,b);
                u2b=du2(:,b);
                u1ab=vec(d2u1(a,b,:));
                u2ab=vec(d2u2(a,b,:));
    
                for c=a:numDerivVars
                    u1c=du1(:,c);
                    u2c=du2(:,c);
                    u1ac=vec(d2u1(a,c,:));
                    u1bc=vec(d2u1(b,c,:));
                    u2ac=vec(d2u2(a,c,:));
                    u2bc=vec(d2u2(b,c,:));
                    u1abc=vec(d3u1(a,b,c,:));
                    u2abc=vec(d3u2(a,b,c,:));
    
                    M(a,b,c)=dot(u1abc,u2)+dot(u1ab,u2c)+dot(u1ac,u2b)+dot(u1a,u2bc)+dot(u1bc,u2a)+dot(u1b,u2ac)+dot(u1c,u2ab)+dot(u1,u2abc);
    
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
