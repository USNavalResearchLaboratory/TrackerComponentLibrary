function [rho,Psi,w]=MEstFunction(u,funcType,a,b,c)
%%MESTFUNC M estimators are maximum-likelihood-like estimators that are
%          typically used for robust linear or nonlinear estimation. Unlike
%          directly performing maximum likelihood esitmation, M estimators
%          transform an error residual term. This function implements a
%          number of common transformation functions used in M estimators
%          and their first derivative as well as a transformed version fo
%          the first derivative that arises in least squares estimation
%          problems.
%
%INPUTS: u A matrix or hypermatrix of real points at which the M estimator
%          value is desired.
% funcType This selects the function to use. The functions listed here are
%          from Table 25.1 of [1]. All functions except 0 use the a input,
%          and function 5 also uses the b and c inputs. Default values for
%          those inputs (if a,b,c are omitted or empty matrices are passed)
%          are shown. Possible values are:
%          0 The least squares criterion.
%          1 (The default if omitted or an empty matrix is passed) Huber's
%            function. a=2, by default.
%          2 Ramsay's function. a=0.3 by default.
%          3 Andrew's wave function with breakpoint a*pi. a=1.339 by
%            default.
%          4 Tukey's biweight function. a=5 by default (suggested in [1] is
%            5<=a<=6).
%          5 Use Hampel's function with breakpoints a, b, and c. a=1.7,
%            b=3.4, and c=8.5 by default.
%          6 Get weights only (rho and Psi are set to empty matrices) based
%            on a weighting function w that is given in Table 25.1 of [1]
%            and is inspired by the t distribution with a degrees of
%            freedom. w(u)=(a+1)/(a+u^2). a=2 by default.
%    a,b,c The breakpoint values used in the selected function. Default
%          values are given above if these are omitted or empty matrices
%          are passed.
%
%OUTPUTS: rho The valus of the selected function at the points in u. 
%         Psi The valus of the derivative of the selected function at the
%             points in u. 
%           w This is Psi./abs(u). These weights arise in linear least
%             squares M estimation. 
%
%Chapter 25.2 of [1] discusses M estimation and these function, which are
%given in Table 25.1.
%
%EXAMPLE:
%This reproduces the plots of Figure 25.1 of [1], except it puts all the
%lines on the same graph. It is plotting w with the default parameters for
%all of the functions.
% numPoints=10000;
% u=linspace(0,6,numPoints).';
% figure(1)
% clf
% hold on
% for k=0:5
%     [~,~,w]=MEstFunction(u,k);
%     plot(u,w,'linewidth',2)
% end
% legend('0','1','2','3','4','5')
% h1=xlabel('u');
% h2=ylabel('w(u)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] N. R. Draper and H. Smith, Applied Regression Analysis, 3rd ed. New
%    York: John Wiley and Sons, Inc., 1998.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(funcType))
    funcType=1;%Huber
end

if(nargin<3||isempty(a))
    switch(funcType)
        case 0%Least squares
        case 1%Huber
            a=2;
        case 2%Ramsay
            a=0.3;
        case 3%Andrew's wave
            a=1.339;
        case 4%Tukey's biweight
            a=5;
        case 5%Hampel
            a=1.7;
        case 6
            a=2;
        otherwise
        error('Unknown function type specified.')
    end
end

%The other two parameters
if(nargin<4||isempty(b))
    b=3.4;
end
if(nargin<5||isempty(c))
    c=8.5;
end

uDims=size(u);
numU=numel(u);
switch(funcType)
    case 0%Least squares
        rho=(1/2)*u.^2;
        Psi=u;
        w=ones(uDims);
    case 1%Huber
        rho=zeros(uDims);
        Psi=zeros(uDims);
        w=zeros(uDims);
        
        for k=1:numU
            uCur=u(k);
            absU=abs(uCur);
            if(absU<=a)
                rho(k)=(1/2)*uCur.^2;
                Psi(k)=uCur;
                w(k)=1;
            else
                rho(k)=a*absU-(1/2)*a^2;
                Psi(k)=a*sign(uCur);
                w(k)=a/absU;
            end
        end
    case 2%Ramsay
        aAbsU=abs(u);
        expVal=exp(-aAbsU);
        
        rho=(1-expVal.*(1+aAbsU))/a^2;
        Psi=u.*expVal;
        w=expVal;
    case 3%Andrew's wave
        rho=zeros(uDims);
        Psi=zeros(uDims);
        w=zeros(uDims);
        
        aPi=a*pi;
        for k=1:numU
            uCur=u(k);
            absU=abs(uCur);
            if(absU<=aPi)
                uda=uCur/a;
                
                rho(k)=a*(1-cos(uda));
                Psi(k)=sin(uda);
                w(k)=Psi(k)/uda;
            else
                rho(k)=2*a;
                Psi(k)=0;
                w(k)=0;
            end
        end
    case 4%Tukey's biweight
        rho=zeros(uDims);
        Psi=zeros(uDims);
        w=zeros(uDims);
        
        for k=1:numU
            uCur=u(k);
            absU=abs(u(k));
            
            if(absU<=a)
                rho(k)=(1/2)*uCur^2-uCur^4/(4*a^2);
                w(k)=1-uCur^2/a^2;
                
                Psi(k)=uCur*w(k);
            else
                rho(k)=(1/4)*a^2;
                Psi(k)=0;
                w(k)=0;
            end
        end
    case 5%Hampel
        rho=zeros(uDims);
        Psi=zeros(uDims);
        w=zeros(uDims);
        
        for k=1:numU
            uCur=u(k);
            absU=abs(u(k));
            
            if(absU<=a)
                rho(k)=(1/2)*uCur^2;
                Psi(k)=uCur;
                w(k)=1;
            elseif(absU<=b)
                rho(k)=a*absU-(1/2)*a^2;
                Psi(k)=a*sign(uCur);
                w(k)=a/absU;
            elseif(absU<=c)
                rho(k)=a*(c*absU-(1/2)*uCur^2)/(c-b)-(7/6)*a^2;
                Psi(k)=(c*sign(uCur)-uCur)*a/(c-b);
                w(k)=(c/absU-1)*a/(c-b);
            else
                rho(k)=a*(b+c-a);
                Psi(k)=0;
                w(k)=0;
            end
        end
    case 6
        rho=[];
        Psi=[];
        w=(a+1)./(a+u.^2);
    otherwise
        error('Unknown function type specified.')
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
