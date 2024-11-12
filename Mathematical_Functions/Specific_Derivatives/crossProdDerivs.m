function [J,H,M]=crossProdDerivs(u1,u2,du1,du2,d2u1,d2u2,d3u1,d3u2)
%%CROSSPRODDERIVS Assuming that the 3X1 vectors u1 and u2 are functions of
%   some other variables, this function can evaluate the first through
%   and up to the third partial derivatives of cross(u1,u2) given partial
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
%OUTPUTS: J A 3XN set of first derivatives of cross(u1,u2), with J(:,i)
%           holding the derivative with respect to the ith parameter.
%         H An NXNX3 set of second derivatives of cross(u1,u2) with
%           H(i,j,:) holding the second derivative with respect to
%           parameters i and j.
%         M An NXNXNX3 set of third derivatives of cross(u1,u2) with
%           M(i,j,k,:) holding the third derivatives with respect to
%           parameters i, j, and k.
%
%This function just returns analytic first second and third deriavtives of
%cross(u,1,u2), where u1 and u2 depend on some variables.
%
%EXAMPLE 1:
%In this example, the first derivatives returned by this function are
%compared to finite differencing for the case where the two vectors are
%unit vectors given in spherical cooridnates. The results agree to more
%than 7 digits, which is good.
% angles=deg2rad([25,-45;
%                 30,10]);
% epsVal=1e-8;
% u=zeros(3,2);
% du=zeros(3,4,2);
% for curVec=1:2
%     a=angles(1,curVec);
%     b=angles(2,curVec);
%     sinA=sin(a);
%     sinB=sin(b);
%     cosA=cos(a);
%     cosB=cos(b);
%     u(:,curVec)=[cosA*cosB;
%                  sinA*cosB;
%                       sinB];
%     %First derivatives
%     idx=(curVec-1)*2;
%     du(:,idx+1,curVec)=[-cosB*sinA;
%            cosA*cosB;
%                    0];
%     du(:,idx+2,curVec)=[-cosA*sinB;
%           -sinA*sinB;
%                 cosB];
% end
% 
% uCross=cross(u(:,1),u(:,2));
% uCrossEps=zeros(3,4);
% for curEpsVar=1:4
%     anglesCur=angles;
%     anglesCur(curEpsVar)=anglesCur(curEpsVar)+epsVal;
% 
%     uVecs=zeros(3,2);
%     for curVec=1:2
%         a=anglesCur(1,curVec);
%         b=anglesCur(2,curVec);
%         sinA=sin(a);
%         sinB=sin(b);
%         cosA=cos(a);
%         cosB=cos(b);
% 
%         uVecs(:,curVec)=[cosA*cosB;
%                          sinA*cosB;
%                               sinB];
%     end
%     uCrossEps(:,curEpsVar)=cross(uVecs(:,1),uVecs(:,2));
% end
% J=crossProdDerivs(u(:,1),u(:,2),du(:,:,1),du(:,:,2));
% JNumDiff=zeros(3,4);
% for k=1:4
%     JNumDiff(:,k)=(uCrossEps(:,k)-uCross)/epsVal;
% end
% RelError=max(abs((JNumDiff(:)-J(:))./JNumDiff(:)))
%
%EXAMPLE 2:
%In this example, the second derivatives returned by this function are
%compared to finite differencing for the case where the two vectors are
%unit vectors given in spherical cooridnates. The results agree to more
%than 7 digits, which is good.
% angles=deg2rad([25,-45;
%                 30,10]);
% epsVal=1e-8;
% u=zeros(3,2);
% du=zeros(3,4,2);
% d2u=zeros(4,4,3,2);
% for curVec=1:2
%     a=angles(1,curVec);
%     b=angles(2,curVec);
%     sinA=sin(a);
%     sinB=sin(b);
%     cosA=cos(a);
%     cosB=cos(b);
%     u(:,curVec)=[cosA*cosB;
%                  sinA*cosB;
%                       sinB];
%     %First derivatives
%     idx=(curVec-1)*2;
%     du(:,idx+1,curVec)=[-cosB*sinA;
%                          cosA*cosB;
%                                  0];
%     du(:,idx+2,curVec)=[-cosA*sinB;
%                         -sinA*sinB;
%                               cosB];
%     %Second derivatives
%     d2u(idx+1,idx+1,:,curVec)=[-cosA*cosB;
%                                -cosB*sinA;
%                                         0];
%     d2u(idx+1,idx+2,:,curVec)=[sinA*sinB;
%                               -cosA*sinB;
%                                        0];
%     d2u(idx+2,idx+1,:,curVec)=d2u(idx+1,idx+2,:,curVec);
%     d2u(idx+2,idx+2,:,curVec)=[-cosA*cosB;
%                                -cosB*sinA;
%                                     -sinB];
% end
% 
% [J,H]=crossProdDerivs(u(:,1),u(:,2),du(:,:,1),du(:,:,2),d2u(:,:,:,1),d2u(:,:,:,2));
% JEps=zeros(3,4,4);
% for curEpsVar=1:4
%     anglesCur=angles;
%     anglesCur(curEpsVar)=anglesCur(curEpsVar)+epsVal;
%     uEps=zeros(3,2);
%     duEps=zeros(3,4,2);
% 
%     for curVec=1:2
%         a=anglesCur(1,curVec);
%         b=anglesCur(2,curVec);
%         sinA=sin(a);
%         sinB=sin(b);
%         cosA=cos(a);
%         cosB=cos(b);
%         uEps(:,curVec)=[cosA*cosB;
%                      sinA*cosB;
%                           sinB];
%         %First derivatives
%         idx=(curVec-1)*2;
%         duEps(:,idx+1,curVec)=[-cosB*sinA;
%                                 cosA*cosB;
%                                         0];
%         duEps(:,idx+2,curVec)=[-cosA*sinB;
%                                -sinA*sinB;
%                                     cosB];
%     end
%     JEps(:,:,curEpsVar)=crossProdDerivs(uEps(:,1),uEps(:,2),duEps(:,:,1),duEps(:,:,2));
% end
% 
% HNumDiff=zeros(4,4,3);
% for k1=1:4
%     JCur=J(:,k1);
%     for k2=1:4
%         JEpsCur=JEps(:,k1,k2);
%         HNumDiff(k1,k2,:)=(JEpsCur-JCur)/epsVal;
%     end
% end
% RelError=max(abs((HNumDiff(:)-H(:))./HNumDiff(:)))
%
%EXAMPLE 3:
%In this example, the third derivatives returned by this function are
%compared to finite differencing for the case where the two vectors are
%unit vectors given in spherical cooridnates. The results agree to more
%than 7 digits, which is good.
% angles=deg2rad([25,-45;
%                 30,10]);
% epsVal=1e-8;
% u=zeros(3,2);
% du=zeros(3,4,2);
% d2u=zeros(4,4,3,2);
% d3u=zeros(4,4,4,3,2);
% for curVec=1:2
%     a=angles(1,curVec);
%     b=angles(2,curVec);
%     sinA=sin(a);
%     sinB=sin(b);
%     cosA=cos(a);
%     cosB=cos(b);
%     u(:,curVec)=[cosA*cosB;
%                  sinA*cosB;
%                       sinB];
%     %First derivatives
%     idx=(curVec-1)*2;
%     du(:,idx+1,curVec)=[-cosB*sinA;
%                          cosA*cosB;
%                                  0];
%     du(:,idx+2,curVec)=[-cosA*sinB;
%                         -sinA*sinB;
%                               cosB];
%     %Second derivatives
%     d2u(idx+1,idx+1,:,curVec)=[-cosA*cosB;
%                                -cosB*sinA;
%                                         0];
%     d2u(idx+1,idx+2,:,curVec)=[sinA*sinB;
%                               -cosA*sinB;
%                                        0];
%     d2u(idx+2,idx+1,:,curVec)=d2u(idx+1,idx+2,:,curVec);
%     d2u(idx+2,idx+2,:,curVec)=[-cosA*cosB;
%                                -cosB*sinA;
%                                     -sinB];
%     %Third derivatives
%     d3u(idx+1,idx+1,idx+1,:,curVec)=[cosB*sinA;
%                                     -cosA*cosB;
%                                              0];
%     d3u(idx+1,idx+1,idx+2,:,curVec)=[cosA*sinB;
%                                      sinA*sinB;
%                                              0];
%     d3u(idx+1,idx+2,idx+2,:,curVec)=[cosB*sinA;
%                                     -cosA*cosB;
%                                              0];
%     d3u(idx+2,idx+2,idx+2,:,curVec)=[cosA*sinB;
%                                      sinA*sinB;
%                                          -cosB];
% 
%     %All of the permutations of the indices are the same.
%     d3u(idx+1,idx+1,idx+2,:,curVec)=d3u(idx+1,idx+1,idx+2,:,curVec);
%     d3u(idx+1,idx+2,idx+1,:,curVec)=d3u(idx+1,idx+1,idx+2,:,curVec);
%     d3u(idx+2,idx+1,idx+2,:,curVec)=d3u(idx+1,idx+1,idx+2,:,curVec);
%     d3u(idx+1,idx+2,idx+2,:,curVec)=d3u(idx+1,idx+2,idx+2,:,curVec);
%     d3u(idx+2,idx+1,idx+2,:,curVec)=d3u(idx+1,idx+2,idx+2,:,curVec);
%     d3u(idx+2,idx+2,idx+1,:,curVec)=d3u(idx+1,idx+2,idx+2,:,curVec);
% end
% [J,H,M]=crossProdDerivs(u(:,1),u(:,2),du(:,:,1),du(:,:,2),d2u(:,:,:,1),d2u(:,:,:,2),d3u(:,:,:,:,1),d3u(:,:,:,:,2));
% 
% HEps=zeros(4,4,3,4);
% for curEpsVar=1:4
%     anglesCur=angles;
%     anglesCur(curEpsVar)=anglesCur(curEpsVar)+epsVal;
%     uEps=zeros(3,2);
%     duEps=zeros(3,4,2);
%     d2uEps=zeros(4,4,3,2);
%     for curVec=1:2
%         a=anglesCur(1,curVec);
%         b=anglesCur(2,curVec);
%         sinA=sin(a);
%         sinB=sin(b);
%         cosA=cos(a);
%         cosB=cos(b);
%         uEps(:,curVec)=[cosA*cosB;
%                      sinA*cosB;
%                           sinB];
%         %First derivatives
%         idx=(curVec-1)*2;
%         duEps(:,idx+1,curVec)=[-cosB*sinA;
%                                 cosA*cosB;
%                                         0];
%         duEps(:,idx+2,curVec)=[-cosA*sinB;
%                                -sinA*sinB;
%                                     cosB];
% 
%         %Second derivatives
%         d2uEps(idx+1,idx+1,:,curVec)=[-cosA*cosB;
%                                    -cosB*sinA;
%                                             0];
%         d2uEps(idx+1,idx+2,:,curVec)=[sinA*sinB;
%                                   -cosA*sinB;
%                                            0];
%         d2uEps(idx+2,idx+1,:,curVec)=d2uEps(idx+1,idx+2,:,curVec);
%         d2uEps(idx+2,idx+2,:,curVec)=[-cosA*cosB;
%                                    -cosB*sinA;
%                                         -sinB];
%     end
%     [~,HEps(:,:,:,curEpsVar)]=crossProdDerivs(uEps(:,1),uEps(:,2),duEps(:,:,1),duEps(:,:,2),d2uEps(:,:,:,1),d2uEps(:,:,:,2));
% end
% 
% MNumDiff=zeros(4,4,4,3);
% for k1=1:4
%     for k2=1:4
%         HCur=H(k1,k2,:);
%         for k3=1:4
%             HEpsCur=HEps(k1,k2,:,k3);
%             MNumDiff(k1,k2,k3,:)=(HEpsCur-HCur)/epsVal;
%         end
%     end
% end
% RelError=max(abs((MNumDiff(:)-M(:))./MNumDiff(:)))
%
%July 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(u1,1);

if(xDim~=3)
    error('Only 3D cross products are supported.')
end

numDerivVars=size(du1,2);

J=zeros(3,numDerivVars);
for a=1:numDerivVars
    u1a=du1(:,a);
    u2a=du2(:,a);

    J(:,a)=cross(u1a,u2)+cross(u1,u2a);
end

if(nargout>1)
    H=zeros(numDerivVars,numDerivVars,3);
    for a=1:numDerivVars
        u1a=du1(:,a);
        u2a=du2(:,a);
    
        for b=a:numDerivVars
            u1b=du1(:,b);
            u2b=du2(:,b);
            u1ab=vec(d2u1(a,b,:));
            u2ab=vec(d2u2(a,b,:));
    
            H(a,b,:)=reshape(cross(u1ab,u2)+cross(u1a,u2b)+cross(u1b,u2a)+cross(u1,u2ab),[1,1,3]);
    
            %Using symmetry.
            H(b,a,:)=H(a,b,:);
        end
    end

    if(nargout>2)
        M=zeros(numDerivVars,numDerivVars,numDerivVars,3);
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
    
                    M(a,b,c,:)=reshape(cross(u1abc,u2)+cross(u1ab,u2c)+cross(u1ac,u2b)+cross(u1a,u2bc)+cross(u1bc,u2a)+cross(u1b,u2ac)+cross(u1c,u2ab)+cross(u1,u2abc),[1,1,1,3]);
    
                    %All of the permutations of the indices are the same.
                    M(b,a,c,:)=M(a,b,c,:);
                    M(c,a,b,:)=M(a,b,c,:);
                    M(a,c,b,:)=M(a,b,c,:);
                    M(b,c,a,:)=M(a,b,c,:);
                    M(c,b,a,:)=M(a,b,c,:);
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
