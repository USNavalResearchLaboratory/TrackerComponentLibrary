function xEst=DopplerOnlyInit6D(rDot,lRx,algorithm,opts,AbsTol,scratchFolderPath,execPath)
%%DOPPLERONLYINIT6D Determine the 3D position and velocity of a target
%          using six simultaneous Doppler measurements. The measurements
%          are assumed to be one-way (from a moving emitter). The receivers
%          are stationary. This function makes use of the function
%          solvePolySysWithExtProg to call an external program to solve the
%          simultaneous multivariate polynomials that define the solution.
%
%INPUTS: rDot A 6X1 or 1X6 vector of Doppler measurements of the emitter.
%     lRx A 3X6 set of the six Cartesian locations of the receivers
%         taking the Doppler measurements in rDot.
% algorithm An optional parameter specifying which algorithm to use. The
%         solvers for simultaneous multivariate polynomials are external
%         programs called via solvePolySysWithExtProg. Possible values
%         are:
%         0 (The default if omitted or an empty matrix is passed) Use
%           Bertini.
%         1 Use PHCpack.
%         2 Use the certified homotopy algorithm that is built into
%           Macaulay2 (in NAG4M2). This only uses the normalized solver
%           with the default options.
%    opts An optional input specifying options for the solver in the
%         function solvePolySysWithExtProg. This is described in more
%         detail in solvePolySysWithExtProg. Omitting this parameter or
%         passing am empty matrices uses the default values, except if
%         Bertini is used, then SecurityMaxNorm, EndpointFiniteThreshold
%         and PathTruncationThreshold are set to 1e9 by default.
%  AbsTol An absolute tolerance on the imaginary part of a solution to the
%         multivariate polynomials used here real. The default if this
%         parameter is omitted or an empty matrix is passed is 1e-8.
% scratchFolderPath An optional parameter specifying the folder in which
%         temporary files will be written. This is needed, because all
%         solvers only works through files. If this parameter is omitted or
%         an empty matrix is passed, then a folder named temp in the folder
%         enclosing this function is used. Note that if an error occurs
%         while executing this function, then temporary files might be left
%         in the temp folder after this function ends.
%  execPath The command line command to use to execute the solver. The 
%         default if this parameter is omitted or an empty matrix is 
%         passed is just the standard name on the command line: bertini,
%         phc and M2 for each of the algorithms. The default assumes that
%         the program is already in the default search path. For example,
%         on *NIX systems, one can usually put the executable in
%         /usr/local/bin for it to be in the default search path.
%         However, for some reason, that directory is often not
%         included in Matlab's search path. Matlab 2016a's help under
%         "Run External Commands, Scripts, and Programs" specifically
%         says how to add that folder to the default search path.
%
%OUTPUTS: xEst A 6XnumEst set of solutions.
%
%The algorithm solves a system of simultaneous multivariate polynomials as
%formulated in [1] (Equation numbers in the code refer to [1]). However,
%there is a typo in [1], which is corrected in the same derivation in [2].
%
%EXAMPLE:
% lRx=[1000,  500, 1100, 2500,     0, 1000;
%      3000, 2500, 2500,    0,     0, 1000;
%         0,    0,    0,  400,  8000, 8000];%(Stationary) sensor locations
% xTrue=[1e3;5e3;4e3;0;300;0];
% rDot=zeros(6,1);
% for curMeas=1:6
%     diff=xTrue(1:3)-lRx(:,curMeas);
%     rDot(curMeas)=(diff'/norm(diff))*xTrue(4:6);
% end
% DopplerOnlyInit6D(rDot,lRx,0)
%
%REFERENCES:
%[1] T.-L. Lee, S.-S. Lin, W.-W. Lin, S.-T. Yau, and J. Zhu, "Polynomial
%    calculations in Doppler tracking," Communications in Information and
%    Systems, vol. 12, no. 2, pp. 157-184, 2012.
%[2] Y.-C. Kuo, W.-W. Lin, and S.-T. Yau, "A novel efficient homotopy
%    continuation method in tracking," Communications in Information and
%    Systems, vol. 14, no. 1, pp. 57-78, 2014.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%We subtract out the locations of the sensors and then put them back in at
%the end. This should help with numerical sensitivity issues when
%localizing things on the ground but using ECEF coordinates. The offset is
%removed in the end.
centerOfRegion=mean(lRx,2);
lRx=bsxfun(@minus,lRx,centerOfRegion);

%Now, we scale the distances so that the farthest one is 1. This should
%help deal with convergence problems when the scale of the probolem
%changes. The scaling will be reversed after the estimation.
scalFactor=sqrt(max(sum(lRx.*lRx,1)));
lRx=lRx/scalFactor;
rDot=rDot/scalFactor;

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

if(nargin<4)
    opts=[];
end

if(nargin<5||isempty(AbsTol))
    AbsTol=1e-8;
end 

if(nargin<5)
    scratchFolderPath=[];
end

if(nargin<6)
    execPath=[];
end

if(algorithm==0&&isempty(opts))%If using Bertini, set default tolerances
    %Default options when using Bertini
    opts=struct('SecurityMaxNorm',1e9,'EndpointFiniteThreshold',1e9,'PathTruncationThreshold',1e9);
elseif(algorithm==0)
    if(~isfield(opts,'SecurityMaxNorm'))
       opts.SecurityMaxNorm=1e9;
    end
    
    if(~isfield(opts,'EndpointFiniteThreshold'))
        opts.EndpointFiniteThreshold=1e9;
    end
    
    if(~isfield(opts,'PathTruncationThreshold'))
        opts.PathTruncationThreshold=1e9;
    end
end

%The following are from Equation 4.
V0=lRx';
n=sum(V0.*V0,2);

u1=lRx(:,1);
u2=lRx(:,2);
u3=lRx(:,3);
u4=lRx(:,4);
u5=lRx(:,5);
u6=lRx(:,6);

%From Equation 5
RDot=diag(rDot);

%From Equation 6
C=[-1, 1, 0, 0, 0, 0;
   -1, 0, 1, 0, 0, 0;
   -1, 0, 0, 1, 0, 0;
   -1, 0, 0, 0, 1, 0;
   -1, 0, 0, 0, 0, 1];

%From Equation 10
A=(V0'*(C'*C)*V0)\(V0'*(C'*C));

%From Equation 15
CHat=[0, -1, 1, 0, 0, 0;
      0, -1, 0, 1, 0, 0;
      0, -1, 0, 0, 1, 0;
      0, -1, 0, 0, 0, 1];

%From Equation 16 --Note that the RDot matrix has been removed. This is
%because including RDot makes the solution inconsistent with Equation 22.
%However, the RDot term needs to be included in the QR decomposition.
AHat=([(u2-u3)';
      (u2-u4)';
      (u2-u5)';
      (u2-u6)']*A+CHat);
  
%Equation 18
u1Hat=(1/2)*A*n-u1;

%Get the T values used in Equation 17 (needs the RDot term with AHat.
[~,R]=qr(AHat*RDot);
T1=R(1:2,1:2);
T2=R(1:2,3:end);

%This is from Equation 17
r12Transform=-T1\T2;
%Extract the elements of the transformation
at11=r12Transform(1,1);
at12=r12Transform(1,2);
at13=r12Transform(1,3);
at14=r12Transform(1,4);
at21=r12Transform(2,1);
at22=r12Transform(2,2);
at23=r12Transform(2,3);
at24=r12Transform(2,4);

%A term in the first formula in Equation 23.
M=(1/2)*(A'*A)*RDot;
%Extract the elements
m1=M(1);
m2=M(2);
m3=M(3);
m4=M(4);
m5=M(5);
m6=M(6);
m7=M(7);
m8=M(8);
m9=M(9);
m10=M(10);
m11=M(11);
m12=M(12);
m13=M(13);
m14=M(14);
m15=M(15);
m16=M(16);
m17=M(17);
m18=M(18);
m19=M(19);
m20=M(20);
m21=M(21);
m22=M(22);
m23=M(23);
m24=M(24);
m25=M(25);
m26=M(26);
m27=M(27);
m28=M(28);
m29=M(29);
m30=M(30);
m31=M(31);
m32=M(32);
m33=M(33);
m34=M(34);
m35=M(35);
m36=M(36);

%A vector term in the first formula in Equation 23.
mVec=u1Hat'*A*RDot;
%Extract the elements
mv1=mVec(1);
mv2=mVec(2);
mv3=mVec(3);
mv4=mVec(4);
mv5=mVec(5);
mv6=mVec(6);

%The first polynomial in Equation 23. The components are ordered
%[r3;r4;r5;r6].
terms1=zeros(5,24);%Allocate space for the terms.

terms1(1,1)=(-at11*(rDot(1)+mv1)-at21*mv2-mv3);
terms1(2:end,1)=[1;0;0;0];%r3

terms1(1,2)=(at11^3*m1+m15+at11*(at21^2*m2+m3)+at11^2*(m13+at21*m7)+at21*(at21*(m14+at21*m8)+m9));
terms1(2:end,2)=[3;0;0;0];%r3^3

terms1(1,3)=(-at12*(rDot(1)+mv1)-at22*mv2-mv4);
terms1(2:end,3)=[0;1;0;0];%r4

terms1(1,4)=(2*at21*at22*m14+m21+at12*m3+at11^2*(3*at12*m1+m19+at22*m7)+2*at11*(at21*at22*m2+at12*(m13+at21*m7))+at21^2*(at12*m2+m20+3*at22*m8)+at22*m9);
terms1(2:end,4)=[2;1;0;0];%r3^2*r4

terms1(1,5)=(at12^2*m13+at22^2*m14+m16+at11*(3*at12^2*m1+at22^2*m2+m4+2*at12*(m19+at22*m7))+at21*(m10+2*at12*at22*m2+2*at22*m20+at12^2*m7+3*at22^2*m8));
terms1(2:end,5)=[1;2;0;0];%r3*r4^2

terms1(1,6)=(at12^3*m1+m22+at12*(at22^2*m2+m4)+at12^2*(m19+at22*m7)+at22*(m10+at22*(m20+at22*m8)));
terms1(2:end,6)=[0;3;0;0];%r4^3

terms1(1,7)=(-at13*(rDot(1)+mv1)-at23*mv2-mv5);
terms1(2:end,7)=[0;0;1;0];%r5

terms1(1,8)=(2*at21*at23*m14+m27+at13*m3+at11^2*(3*at13*m1+m25+at23*m7)+2*at11*(at21*at23*m2+at13*(m13+at21*m7))+at21^2*(at13*m2+m26+3*at23*m8)+at23*m9);
terms1(2:end,8)=[2;0;1;0];%r3^2*r5

terms1(1,9)=2*(at22*at23*m14+at12*(at21*at23*m2+at13*(m13+at21*m7))+at11*(at13*m19+at22*at23*m2+at13*at22*m7+at12*(3*at13*m1+m25+at23*m7))+at21*(at23*m20+at22*(at13*m2+m26+3*at23*m8)));
terms1(2:end,9)=[1;1;1;0];%r3*r4*r5

terms1(1,10)=(at22^2*(at13*m2+m26)+m28+at13*m4+at12^2*(3*at13*m1+m25+at23*m7)+2*at12*(at22*at23*m2+at13*(m19+at22*m7))+at23*(m10+2*at22*m20+3*at22^2*m8));
terms1(2:end,10)=[0;2;1;0];%r4^2*r5

terms1(1,11)=(at13^2*m13+at23^2*m14+m17+at11*(3*at13^2*m1+at23^2*m2+m5+2*at13*(m25+at23*m7))+at21*(m11+2*at13*at23*m2+2*at23*m26+at13^2*m7+3*at23^2*m8));
terms1(2:end,11)=[1;0;2;0];%r3*r5^2

terms1(1,12)=(at22*m11+2*at13*at22*at23*m2+m23+at13^2*(m19+at22*m7)+at12*(3*at13^2*m1+at23^2*m2+m5+2*at13*(m25+at23*m7))+at23*(2*at22*m26+at23*(m20+3*at22*m8)));
terms1(2:end,12)=[0;1;2;0];%r4*r5^2;

terms1(1,13)=(at13^3*m1+m29+at13*(at23^2*m2+m5)+at13^2*(m25+at23*m7)+at23*(m11+at23*(m26+at23*m8)));
terms1(2:end,13)=[0;0;3;0];%r5^3

terms1(1,14)=(-at14*(rDot(1)+mv1)-at24*mv2-mv6);
terms1(2:end,14)=[0;0;0;1];%r6

terms1(1,15)=(2*at21*at24*m14+at14*m3+m33+at11^2*(3*at14*m1+m31+at24*m7)+2*at11*(at21*at24*m2+at14*(m13+at21*m7))+at21^2*(at14*m2+m32+3*at24*m8)+at24*m9);
terms1(2:end,15)=[2;0;0;1];%r3^2*r6

terms1(1,16)=2*(at22*at24*m14+at12*(at21*at24*m2+at14*(m13+at21*m7))+at11*(at14*m19+at22*at24*m2+at14*at22*m7+at12*(3*at14*m1+m31+at24*m7))+at21*(at24*m20+at22*(at14*m2+m32+3*at24*m8)));
terms1(2:end,16)=[1;1;0;1];%r3*r4*r6

terms1(1,17)=(at22^2*(at14*m2+m32)+m34+at14*m4+at12^2*(3*at14*m1+m31+at24*m7)+2*at12*(at22*at24*m2+at14*(m19+at22*m7))+at24*(m10+2*at22*m20+3*at22^2*m8));
terms1(2:end,17)=[0;2;0;1];%r4^2*r6

terms1(1,18)=2*(at23*at24*m14+at13*(at21*at24*m2+at14*(m13+at21*m7))+at11*(at23*at24*m2+at14*m25+at14*at23*m7+at13*(3*at14*m1+m31+at24*m7))+at21*(at24*m26+at23*(at14*m2+m32+3*at24*m8)));
terms1(2:end,18)=[1;0;1;1];%r3*r5*r6

terms1(1,19)=2*(at14*at22*at23*m2+at23*at24*m20+at22*at24*m26+at22*at23*m32+at13*(at14*m19+at22*at24*m2+at14*at22*m7)+at12*(at23*at24*m2+at14*m25+at14*at23*m7+at13*(3*at14*m1+m31+at24*m7))+3*at22*at23*at24*m8);
terms1(2:end,19)=[0;1;1;1];%r4*r5*r6

terms1(1,20)=(at23^2*(at14*m2+m32)+m35+at14*m5+2*at13*(at23*at24*m2+at14*m25+at14*at23*m7)+at13^2*(3*at14*m1+m31+at24*m7)+at24*(m11+2*at23*m26+3*at23^2*m8));
terms1(2:end,20)=[0;0;2;1];%r5^2*r6

terms1(1,21)=(at14^2*m13+at24^2*m14+m18+at11*(3*at14^2*m1+at24^2*m2+m6+2*at14*(m31+at24*m7))+at21*(m12+2*at14*at24*m2+2*at24*m32+at14^2*m7+3*at24^2*m8));
terms1(2:end,21)=[1;0;0;2];%r3*r6^2

terms1(1,22)=(at22*m12+2*at14*at22*at24*m2+m24+at14^2*(m19+at22*m7)+at12*(3*at14^2*m1+at24^2*m2+m6+2*at14*(m31+at24*m7))+at24*(2*at22*m32+at24*(m20+3*at22*m8))) ;
terms1(2:end,22)=[0;1;0;2];%r4*r6^2

terms1(1,23)=(at23*m12+2*at14*at23*at24*m2+m30+at14^2*(m25+at23*m7)+at13*(3*at14^2*m1+at24^2*m2+m6+2*at14*(m31+at24*m7))+at24*(2*at23*m32+at24*(m26+3*at23*m8)));
terms1(2:end,23)=[0;0;1;2];%r5*r6^2

terms1(1,24)=(at14^3*m1+m36+at14*(at24^2*m2+m6)+at14^2*(m31+at24*m7)+at24*(m12+at24*(m32+at24*m8)));
terms1(2:end,24)=[0;0;0;3];%r6^3

%Scale the coefficients
terms1(1,:)=terms1(1,:)/max(abs(terms1(1,:)));

%cHat is defined after Equation 20.
cHat=(1/4)*n'*(A'*A)*n-u1'*A*n+u1'*u1;

%A term in the second formula in Equation 23.
Mt=(1/4)*(A'*A);
%Extract the elements
mt1=Mt(1);
mt2=Mt(2);
mt3=Mt(3);
mt4=Mt(4);
mt5=Mt(5);
mt6=Mt(6);
mt7=Mt(7);
mt8=Mt(8);
mt9=Mt(9);
mt10=Mt(10);
mt11=Mt(11);
mt12=Mt(12);
mt13=Mt(13);
mt14=Mt(14);
mt15=Mt(15);
mt16=Mt(16);
mt17=Mt(17);
mt18=Mt(18);
mt19=Mt(19);
mt20=Mt(20);
mt21=Mt(21);
mt22=Mt(22);
mt23=Mt(23);
mt24=Mt(24);
mt25=Mt(25);
mt26=Mt(26);
mt27=Mt(27);
mt28=Mt(28);
mt29=Mt(29);
mt30=Mt(30);
mt31=Mt(31);
mt32=Mt(32);
mt33=Mt(33);
mt34=Mt(34);
mt35=Mt(35);
mt36=Mt(36);

%A vector term in the second formula in Equation 23.
mVec2=u1Hat'*A;
%Extract the elements
mvt1=mVec2(1);
mvt2=mVec2(2);
mvt3=mVec2(3);
mvt4=mVec2(4);
mvt5=mVec2(5);
mvt6=mVec2(6);

%The second polynomial in Equation 23. The components are ordered
%[r3;r4;r5;r6].
terms2=zeros(5,46);%Allocate space for the terms.

terms2(1,1)=cHat;
terms2(2:end,1)=[0;0;0;0];%Constant term

terms2(1,2)=(at11^4*mt1+mt15+at11^2*(mt13+mt3+at21^2*(mt2+mt7))+at21^4*mt8+at21^2*(mt14+mt9));
terms2(2:end,2)=[4;0;0;0];%r3^4

terms2(1,3)=2*(2*at11^3*at12*mt1+at11^2*at21*at22*(mt2+mt7)+at11*at12*(mt13+mt3+at21^2*(mt2+mt7))+at21*at22*(mt14+2*at21^2*mt8+mt9));
terms2(2:end,3)=[3;1;0;0];%r3^3*r4

terms2(1,4)=(mt16+mt21+at12^2*(mt13+mt3)+4*at11*at12*at21*at22*(mt2+mt7)+at11^2*(6*at12^2*mt1+mt19+mt4+at22^2*(mt2+mt7))+at21^2*(mt10+mt20+at12^2*(mt2+mt7)+6*at22^2*mt8)+at22^2*(mt14+mt9));
terms2(2:end,4)=[2;2;0;0];%r3^2*r4^2

terms2(1,5)=2*(at11*at12*(2*at12^2*mt1+mt19+mt4+at22^2*(mt2+mt7))+at21*at22*(mt10+mt20+at12^2*(mt2+mt7)+2*at22^2*mt8));
terms2(2:end,5)=[1;3;0;0];%r3*r4^3

terms2(1,6)=(at12^4*mt1+at22^2*(mt10+mt20)+mt22+at12^2*(mt19+mt4+at22^2*(mt2+mt7))+at22^4*mt8);
terms2(2:end,6)=[0;4;0;0];%r4^4

terms2(1,7)=2*(2*at11^3*at13*mt1+at11^2*at21*at23*(mt2+mt7)+at11*at13*(mt13+mt3+at21^2*(mt2+mt7))+at21*at23*(mt14+2*at21^2*mt8+mt9));
terms2(2:end,7)=[3;0;1;0];%r3^3*r5

terms2(1,8)=2*(2*at12^3*at13*mt1+at12^2*at22*at23*(mt2+mt7)+at12*at13*(mt19+mt4+at22^2*(mt2+mt7))+at22*at23*(mt10+mt20+2*at22^2*mt8));
terms2(2:end,8)=[0;3;1;0];%r4^3*r5

terms2(1,9)=(mt17+mt27+at13^2*(mt13+mt3)+4*at11*at13*at21*at23*(mt2+mt7)+at11^2*(6*at13^2*mt1+mt25+mt5+at23^2*(mt2+mt7))+at21^2*(mt11+mt26+at13^2*(mt2+mt7)+6*at23^2*mt8)+at23^2*(mt14+mt9));
terms2(2:end,9)=[2;0;2;0];%r3^2*r5^2

terms2(1,10)=(mt23+mt28+at13^2*(mt19+mt4)+4*at12*at13*at22*at23*(mt2+mt7)+at22^2*(mt11+mt26+at13^2*(mt2+mt7))+at12^2*(6*at13^2*mt1+mt25+mt5+at23^2*(mt2+mt7))+at23^2*(mt10+mt20+6*at22^2*mt8));
terms2(2:end,10)=[0;2;2;0];%r4^2*r5^2

terms2(1,11)=2*(at11*at13*(2*at13^2*mt1+mt25+mt5+at23^2*(mt2+mt7))+at21*at23*(mt11+mt26+at13^2*(mt2+mt7)+2*at23^2*mt8));
terms2(2:end,11)=[1;0;3;0];%r3*r5^3

terms2(1,12)=2*(at12*at13*(2*at13^2*mt1+mt25+mt5+at23^2*(mt2+mt7))+at22*at23*(mt11+mt26+at13^2*(mt2+mt7)+2*at23^2*mt8));
terms2(2:end,12)=[0;1;3;0];%r4*r5^3

terms2(1,13)=(at13^4*mt1+at23^2*(mt11+mt26)+mt29+at13^2*(mt25+mt5+at23^2*(mt2+mt7))+at23^4*mt8);
terms2(2:end,13)=[0;0;4;0];%r5^4

terms2(1,14)=-2*(at11*at14*(1+mvt1)+at21*at24*mvt2);
terms2(2:end,14)=[1;0;0;1];%r3*r6

terms2(1,15)=2*(2*at11^3*at14*mt1+at11^2*at21*at24*(mt2+mt7)+at11*at14*(mt13+mt3+at21^2*(mt2+mt7))+at21*at24*(mt14+2*at21^2*mt8+mt9));
terms2(2:end,15)=[3;0;0;1];%r3^3*r6

terms2(1,16)=-2*(at12*at14*(1+mvt1)+at22*at24*mvt2);
terms2(2:end,16)=[0;1;0;1];%r4*r6

terms2(1,17)=2*(2*at12^3*at14*mt1+at12^2*at22*at24*(mt2+mt7)+at12*at14*(mt19+mt4+at22^2*(mt2+mt7))+at22*at24*(mt10+mt20+2*at22^2*mt8));
terms2(2:end,17)=[0;3;0;1];%r4^3*r6

terms2(1,18)=-2*(at13*at14*(1+mvt1)+at23*at24*mvt2);
terms2(2:end,18)=[0;0;1;1];%r5*r6

terms2(1,19)=2*(2*at13^3*at14*mt1+at13^2*at23*at24*(mt2+mt7)+at13*at14*(mt25+mt5+at23^2*(mt2+mt7))+at23*at24*(mt11+mt26+2*at23^2*mt8));
terms2(2:end,19)=[0;0;3;1];%r5^3*r6

terms2(1,20)=(-at14^2*(1+mvt1)-at24^2*mvt2-mvt6);
terms2(2:end,20)=[0;0;0;2];%r6^2

terms2(1,21)=(mt18+at14^2*(mt13+mt3)+mt33+4*at11*at14*at21*at24*(mt2+mt7)+at11^2*(6*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7))+at21^2*(mt12+mt32+at14^2*(mt2+mt7)+6*at24^2*mt8)+at24^2*(mt14+mt9));
terms2(2:end,21)=[2;0;0;2];%r3^2*r6^2

terms2(1,22)=(mt24+mt34+at14^2*(mt19+mt4)+4*at12*at14*at22*at24*(mt2+mt7)+at22^2*(mt12+mt32+at14^2*(mt2+mt7))+at12^2*(6*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7))+at24^2*(mt10+mt20+6*at22^2*mt8));
terms2(2:end,22)=[0;2;0;2];%r4^2*r6^2

terms2(1,23)=(mt30+mt35+at14^2*(mt25+mt5)+4*at13*at14*at23*at24*(mt2+mt7)+at23^2*(mt12+mt32+at14^2*(mt2+mt7))+at13^2*(6*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7))+at24^2*(mt11+mt26+6*at23^2*mt8));
terms2(2:end,23)=[0;0;2;2];%r5^2*r6^2

terms2(1,24)=2*(at11*at14*(2*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7))+at21*at24*(mt12+mt32+at14^2*(mt2+mt7)+2*at24^2*mt8));
terms2(2:end,24)=[1;0;0;3];%r3*r6^3

terms2(1,25)=2*(at12*at14*(2*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7))+at22*at24*(mt12+mt32+at14^2*(mt2+mt7)+2*at24^2*mt8));
terms2(2:end,25)=[0;1;0;3];%r4*r6^3

terms2(1,26)=2*(at13*at14*(2*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7))+at23*at24*(mt12+mt32+at14^2*(mt2+mt7)+2*at24^2*mt8));
terms2(2:end,26)=[0;0;1;3];%r5*r6^3

terms2(1,27)=(at14^4*mt1+at24^2*(mt12+mt32)+mt36+at14^2*(mt31+mt6+at24^2*(mt2+mt7))+at24^4*mt8);
terms2(2:end,27)=[0;0;0;4];%r6^4

terms2(1,28)=(-at13^2*(1+mvt1)-at23^2*mvt2-mvt5);
terms2(2:end,28)=[0;0;2;0];%r5^2

terms2(1,29)=2*(at11*(6*at13^2*at14*mt1+2*at13*at23*at24*(mt2+mt7)+at14*(mt25+mt5+at23^2*(mt2+mt7)))+at21*(2*at13*at14*at23*(mt2+mt7)+at24*(mt11+mt26+at13^2*(mt2+mt7)+6*at23^2*mt8)));
terms2(2:end,29)=[1;0;2;1];%r3*r5^2*r6

terms2(1,30)=2*(at12*(6*at13^2*at14*mt1+2*at13*at23*at24*(mt2+mt7)+at14*(mt25+mt5+at23^2*(mt2+mt7)))+at22*(2*at13*at14*at23*(mt2+mt7)+at24*(mt11+mt26+at13^2*(mt2+mt7)+6*at23^2*mt8)));
terms2(2:end,30)=[0;1;2;1];%r4*r5^2*r6

terms2(1,31)=(-at12^2*(1+mvt1)-at22^2*mvt2-mvt4);
terms2(2:end,31)=[0;2;0;0];%r4^2

terms2(1,32)=2*(at11*(6*at12^2*at13*mt1+2*at12*at22*at23*(mt2+mt7)+at13*(mt19+mt4+at22^2*(mt2+mt7)))+at21*(2*at12*at13*at22*(mt2+mt7)+at23*(mt10+mt20+at12^2*(mt2+mt7)+6*at22^2*mt8)));
terms2(2:end,32)=[1;2;1;0];%r3*r4^2*r5

terms2(1,33)=2*(at11*(6*at12^2*at14*mt1+2*at12*at22*at24*(mt2+mt7)+at14*(mt19+mt4+at22^2*(mt2+mt7)))+at21*(2*at12*at14*at22*(mt2+mt7)+at24*(mt10+mt20+at12^2*(mt2+mt7)+6*at22^2*mt8)));
terms2(2:end,33)=[1;2;0;1];%r3*r4^2*r6

terms2(1,34)=2*(2*at12*at22*(at14*at23+at13*at24)*(mt2+mt7)+at13*at14*(mt19+mt4+at22^2*(mt2+mt7))+at12^2*(6*at13*at14*mt1+at23*at24*(mt2+mt7))+at23*at24*(mt10+mt20+6*at22^2*mt8));
terms2(2:end,34)=[0;2;1;1];%r4^2*r5*r6

terms2(1,35)=(-at11^2*(1+mvt1)-at21^2*mvt2-mvt3);
terms2(2:end,35)=[2;0;0;0];%r3^2

terms2(1,36)=2*(2*at11*at21*(at13*at22+at12*at23)*(mt2+mt7)+at12*at13*(mt13+mt3+at21^2*(mt2+mt7))+at11^2*(6*at12*at13*mt1+at22*at23*(mt2+mt7))+at22*at23*(mt14+6*at21^2*mt8+mt9));
terms2(2:end,36)=[2;1;1;0];%r3^2*r4*r5

terms2(1,37)=2*(2*at11*at21*(at14*at22+at12*at24)*(mt2+mt7)+at12*at14*(mt13+mt3+at21^2*(mt2+mt7))+at11^2*(6*at12*at14*mt1+at22*at24*(mt2+mt7))+at22*at24*(mt14+6*at21^2*mt8+mt9));
terms2(2:end,37)=[2;1;0;1];%r3^2*r4*r6

terms2(1,38)=2*(2*at11*at21*(at14*at23+at13*at24)*(mt2+mt7)+at13*at14*(mt13+mt3+at21^2*(mt2+mt7))+at11^2*(6*at13*at14*mt1+at23*at24*(mt2+mt7))+at23*at24*(mt14+6*at21^2*mt8+mt9));
terms2(2:end,38)=[2;0;1;1];%r3^2*r5*r6

terms2(1,39)=-2*(at11*at12*(1+mvt1)+at21*at22*mvt2);
terms2(2:end,39)=[1;1;0;0];%r3*r4

terms2(1,40)=2*(at11*(2*at13*at22*at23*(mt2+mt7)+at12*(6*at13^2*mt1+mt25+mt5+at23^2*(mt2+mt7)))+at21*(2*at12*at13*at23*(mt2+mt7)+at22*(mt11+mt26+at13^2*(mt2+mt7)+6*at23^2*mt8)));
terms2(2:end,40)=[1;1;2;0];%r3*r4*r5^2

terms2(1,41)=4*(at21*(at13*at14*at22+at12*at14*at23+at12*at13*at24)*(mt2+mt7)+at11*(at22*(at14*at23+at13*at24)*(mt2+mt7)+at12*(6*at13*at14*mt1+at23*at24*(mt2+mt7)))+6*at21*at22*at23*at24*mt8);
terms2(2:end,41)=[1;1;1;1];%r3*r4*r5*r6

terms2(1,42)=2*(at11*(2*at14*at22*at24*(mt2+mt7)+at12*(6*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7)))+at21*(2*at12*at14*at24*(mt2+mt7)+at22*(mt12+mt32+at14^2*(mt2+mt7)+6*at24^2*mt8)));
terms2(2:end,42)=[1;1;0;2];%r3*r4*r6^2

terms2(1,43)=-2*(at11*at13*(1+mvt1)+at21*at23*mvt2);
terms2(2:end,43)=[1;0;1;0];%r3*r5

terms2(1,44)=2*(at11*(2*at14*at23*at24*(mt2+mt7)+at13*(6*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7)))+at21*(2*at13*at14*at24*(mt2+mt7)+at23*(mt12+mt32+at14^2*(mt2+mt7)+6*at24^2*mt8)));
terms2(2:end,44)=[1;0;1;2];%r3*r5*r6^2

terms2(1,45)=-2*(at12*at13*(1+mvt1)+at22*at23*mvt2);
terms2(2:end,45)=[0;1;1;0];%r4*r5

terms2(1,46)=2*(at12*(2*at14*at23*at24*(mt2+mt7)+at13*(6*at14^2*mt1+mt31+mt6+at24^2*(mt2+mt7)))+at22*(2*at13*at14*at24*(mt2+mt7)+at23*(mt12+mt32+at14^2*(mt2+mt7)+6*at24^2*mt8)));
terms2(2:end,46)=[0;1;1;2];%r4*r5*r6^2

%Scale the coefficients
terms2(1,:)=terms2(1,:)/max(abs(terms2(1,:)));

%The third polynomial in Equation 23:
%A vector term in the third formula in Equation 23.
mVec3=(u2-u3)'*A;

%Extract the elements
mvtt1=mVec3(1);
mvtt2=mVec3(2);
mvtt3=mVec3(3);
mvtt4=mVec3(4);
mvtt5=mVec3(5);
mvtt6=mVec3(6);

terms3=zeros(5,11);

terms3(1,1)=(-u3'*u3+u2'*u2-(u2-u3)'*A*n);
terms3(2:end,1)=[0;0;0;0];%Constant term

terms3(1,2)=1+at11^2*mvtt1+at21^2*(-1+mvtt2)+mvtt3;
terms3(2:end,2)=[2;0;0;0];%r3^2

terms3(1,3)=2*(at11*at12*mvtt1+at21*at22*(-1+mvtt2));
terms3(2:end,3)=[1;1;0;0];%r3*r4

terms3(1,4)=(at12^2*mvtt1+at22^2*(-1+mvtt2)+mvtt4);
terms3(2:end,4)=[0;2;0;0];%r4^2

terms3(1,5)=2*(at11*at13*mvtt1+at21*at23*(-1+mvtt2));
terms3(2:end,5)=[1;0;1;0];%r3*r5

terms3(1,6)=2*(at12*at13*mvtt1+at22*at23*(-1+mvtt2));
terms3(2:end,6)=[0;1;1;0];%r4*r5

terms3(1,7)=(at13^2*mvtt1+at23^2*(-1+mvtt2)+mvtt5);
terms3(2:end,7)=[0;0;2;0];%r5^2

terms3(1,8)=2*(at11*at14*mvtt1+at21*at24*(-1+mvtt2));
terms3(2:end,8)=[1;0;0;1];%r3*r6

terms3(1,9)=2*(at12*at14*mvtt1+at22*at24*(-1+mvtt2));
terms3(2:end,9)=[0;1;0;1];%r4*r6

terms3(1,10)=2*(at13*at14*mvtt1+at23*at24*(-1+mvtt2));
terms3(2:end,10)=[0;0;1;1];%r5*r6

terms3(1,11)=(at14^2*mvtt1+at24^2*(-1+mvtt2)+mvtt6);
terms3(2:end,11)=[0;0;0;2];%r6^2

%Scale the coefficients
terms3(1,:)=terms3(1,:)/max(abs(terms3(1,:)));

%The fourth polynomial in Equation 23:
mVec4=(u2-u4)'*A;
%Extract the elements
mvtt1=mVec4(1);
mvtt2=mVec4(2);
mvtt3=mVec4(3);
mvtt4=mVec4(4);
mvtt5=mVec4(5);
mvtt6=mVec4(6);

terms4=zeros(5,11);

terms4(1,1)=(-u4'*u4+u2'*u2-(u2-u4)'*A*n);
terms4(2:end,1)=[0;0;0;0];%Constant term

terms4(1,2)=(at11^2*mvtt1 + at21^2*(-1 + mvtt2) + mvtt3) ;
terms4(2:end,2)=[2;0;0;0];%r3^2

terms4(1,3)=2*(at11*at12*mvtt1 + at21*at22*(-1 + mvtt2)) ;
terms4(2:end,3)=[1;1;0;0];%r3*r4

terms4(1,4)=(1 + at12^2*mvtt1 + at22^2*(-1 + mvtt2) + mvtt4);
terms4(2:end,4)=[0;2;0;0];%r4^2

terms4(1,5)=2*(at11*at13*mvtt1 + at21*at23*(-1 + mvtt2));
terms4(2:end,5)=[1;0;1;0];%r3*r5

terms4(1,6)=2*(at12*at13*mvtt1 + at22*at23*(-1 + mvtt2));
terms4(2:end,6)=[0;1;1;0];%r4*r5

terms4(1,7)=(at13^2*mvtt1 + at23^2*(-1 + mvtt2) + mvtt5);
terms4(2:end,7)=[0;0;2;0];%r5^2

terms4(1,8)=2*(at11*at14*mvtt1 + at21*at24*(-1 + mvtt2));
terms4(2:end,8)=[1;0;0;1];%r3*r6

terms4(1,9)=2*(at12*at14*mvtt1 + at22*at24*(-1 + mvtt2));
terms4(2:end,9)=[0;1;0;1];%r4*r6

terms4(1,10)=2*(at13*at14*mvtt1 + at23*at24*(-1 + mvtt2));
terms4(2:end,10)=[0;0;1;1];%r5*r6

terms4(1,11)=(at14^2*mvtt1 + at24^2*(-1 + mvtt2) + mvtt6);
terms4(2:end,11)=[0;0;0;2];%r6^2

%Scale the coefficients
terms4(1,:)=terms4(1,:)/max(abs(terms4(1,:)));

%Now, put the term matrices into formats that can be used by the
%multivariate polynomial solvers.
xPolys=cell(4,1);
xPolys{1}=terms1;
xPolys{2}=terms2;
xPolys{3}=terms3;
xPolys{4}=terms4;

varNames={'rc','rd','re','rf'};

for curPoly=1:4
    termMat=xPolys{curPoly};

    numTerms=size(termMat,2);

    thePoly=[];

    for curTerm=1:numTerms
        curCoeff=num2str(termMat(1,curTerm),16);

        curMonomial=[];

        for curVar=1:4
            if(termMat(curVar+1,curTerm)~=0)
                if(~isempty(curMonomial))
                   curMonomial=[curMonomial,'*']; 
                end

                if(termMat(curVar+1,curTerm)==1)
                    curMonomial=[curMonomial,varNames{curVar}];
                else
                    curMonomial=[curMonomial,varNames{curVar},'^',num2str(termMat(curVar+1,curTerm))];
                end
            end
        end

        if(isempty(curMonomial))%If it is a constant term.
            %The constant term (if there is one) should be the first
            %term, so we can just set thePoly to it.
            thePoly=curCoeff;
        elseif(~isempty(thePoly))
            thePoly=[thePoly,'+(',curCoeff,')*',curMonomial];
        else%If there is no constant term.
            thePoly=['(',curCoeff,')*',curMonomial];
        end
    end

    xPolys{curPoly}=thePoly;
end

rEst=solvePolySysWithExtProg(xPolys,varNames,algorithm,opts,scratchFolderPath,execPath);

if(isempty(rEst))
    %The solver failed.
    xEst=[];
    return;
end
%Throw out complex solutions. Solutions are deemed complex if the
%imaginary part exceeds AbsTol in magnitude.
sel=all(abs(imag(rEst))<AbsTol,1);
rEst=real(rEst(:,sel));

%Only positive range values are valid.
sel=all(rEst>0,1);
rEst=rEst(:,sel);

numSol=size(rEst,2);
xEst=zeros(6,numSol);

%Given solutions to the r variables, we need to extract the state.
numAdded=0;
for curSol=1:numSol
    %From Equation 17
    r1r2=r12Transform*rEst(:,curSol);
    
    %Only positive ranges are valid.
    if(all(r1r2)>0)    
        r=[r1r2;rEst(:,curSol)];

        %Equation 8
        uDot=A*(-RDot*r);

        %Equation 9
        u=(1/2)*A*(n-r.*r);

        numAdded=numAdded+1;
        xEst(:,numAdded)=[u;uDot];
    end
end

%Shrink to fit the actual number of solutions added.
xEst=xEst(:,1:numAdded);

%Undo the effects of scaling:
xEst=scalFactor*xEst;

%Add back in the offset that was present due to centering the coordinate
%system around the sensors.
xEst(1:3,:)=bsxfun(@plus,xEst(1:3,:),centerOfRegion);
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
