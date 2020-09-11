function [DVal,totalError,exitCode]=KLDistApproxGaussMix(w1,mu1,P1,w2,mu2,P2,algorithm,param8,w)
%%KLDISTAPPROXGAUSSMIX Approximate the Kullback-Leibler (KL) divergence
%          (also known as the KL distance or relative entropy) between two
%          multivariate Gaussian mixture distributions. Whereas an explicit
%          solution is available for the KL distance between two Gaussian
%          PDFs, no explicit solution exists between a Gaussian PDF and a
%          Gaussian mixture. Thus, various approximations are necessary.
%          As defined in Chapter 8.5 of [4], the relative entropy between
%          two continuous PDFs, f and g, is
%          D(f||g)=integral_{x} f(x) log(f(x)/g(x)) dx
%          The integral is over an appropriate domain of x, in this case
%          the d-dimensional real space for an n-dimensional state. Note
%          that changing the order of f and g changes the result (the
%          operation is not commutative).
%
%INPUTS: w1 The 1Xn1 or n1X1 set of weights of the first Gaussian mixture's
%           components such that all w>=0 and sum(w)=1.
%       mu1 The dXn1 set of mean values of the first Gaussian mixture's
%           components.
%        P1 The dXxdXn1 set of covariance matrices of the first Gaussian
%           mixture's components.
% w2,mu2,P2 The equivalent of w1, mu1 and P1 for the second mixtures
%           components. There length-n2, dXn2, and dXdXn2 in size.
% algorithm An optional parameter specifying the algorithm used to obtain
%           the approximation of the KL divergence. Possible values are
%           for:
%           0 Use adaptive numeric integration to directly evaluate the
%             Kullback-Leiberler divergence integral over a rectangular 
%             region that is standard deviations from the means of the
%             distributions. The default number of standard deviations in
%             each direction is 4.5 and can be adjusted by passing a member
%             of param8 named 'numStdDev'. In 1D, the function
%             integral1DAdaptiveCC is used for the integration and the
%             inputs to that function 'RelTol', 'AbsTol', 'NPowMin', and
%             'NPowMax' can be passed as members of param8 if one wishes to
%             not use the defaults. If higher dimensions, the function
%             integrateUnifCub57Adaptive is used for the numeric
%             integration and the inputs to that function 'maxSearchReg',
%             'AbsTol', and 'RelTol' can be passed if one wishes to not use
%              the defaults. This approximation can produce negative
%              values, though it is often good.
%           1 Monte Carlo sampling from Equation 4 of Section 2 of [3].
%             This uses param8 as the number of samples (if given). The
%             default if param8 is not given is 3e3. As the number of
%             samples increases, this method should become increasingly
%             accurate. This approximation can produce negative values.
%           2 (The default if omitted or an empty matrix is passed) The
%             cubature integration approximation from Equation 8 in Section
%             3 of [3]. This uses param8 as the cubature points and w as
%             the weights, if provided. If not provided, then
%             fifthOrderCubPoints is used to obtain fifth-order xDim-
%             dimensional cubature points. This approximation can produce
%             negative values.
%           3 Use the variational upper bound from Section 8 of [3]. This
%             approximation should always be positive. 10 iterations are
%             performed. This bound is typically not zero when the two
%             input PDFs are equal.
%           4 The variational approximation from Section 7 of [3] (also
%             given in [1]). This approximation can produce negative
%             values.
%           5 The product of Gaussians approximation in [1] and Section 5
%             of [3]. This approximation can produce negative values.
%           6 Take the average of algorithms 5 and 4, which tend to be
%             upper and lower bounds.
%           7 Goldberger's approximation from Equation 2 of [2]. This is
%             also given in Section 6 of [3] as the matched bound
%             approximation. This approximation can produce negative
%             values and is often quite bad with certain classes of
%             distributions.
%    param8 An optional parameter. If algorithm=0, then this is a structure
%           that can take members as described above for algorithm 0. If
%           algorithm=1, then this is the number of samples to use. If
%           algorithm=2, then this is a set of cubature points as a
%           dXnumPoints matrix. For other values of algorithm, this input
%           is not used.
%         w This input is only used if algorithm=2, in which case it is the
%           set of cubature weights (sum(w)=1). that are associated with
%           the cubature points in param8. If param8 is omitted or an empty
%           matrix is passed, then w will be the weights associated with
%           the points from the fifthOrderCubPoints function.
%
%%OUTPUTS: val The approximate value of the Kullback-Leibler divergence
%              between the first and the second Gaussian mixture
%              distributions.
%
%Some approximations can produce negative values (the true KL divergence
%can never be negative).
%
%EXAMPLE 1:
%Though many of the algorithms are often cited in the literature, many of
%the approximations can be quite bad. Here, we consider a 1D example and
%then the mixture missing a few components and also the moment-matched
%Gaussian approximation.
% w1=[0.03,0.18,0.12,0.19,0.02,0.16,0.06,0.1,0.08,0.06];
% n1=length(w1);
% mu1=[1.45,2.20,0.67,0.48,1.49,0.91,1.01,1.42,2.77,0.89];
% P1=[0.0487,0.0305,0.1171,0.0174,0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P1=reshape(P1,[1,1,n1]);
% 
% %The second PDF is the first with the five least-weight components deleted.
% w2=[0.18,0.12,0.19,0.16,0.1,0.08];
% w2=w2/sum(w2);
% n2=length(w2);
% mu2=[2.20,0.67,0.48,0.91,1.42,2.77];
% P2=[0.0305,0.1171,0.0174,0.0102,0.0380,0.0115];
% P2=reshape(P2,[1,1,n2]);
% 
% %The third PDF is the second PDF with two more components deleted.
% w3=[0.18,0.12,0.19,0.16];
% w3=w3/sum(w3);
% n3=length(w3);
% mu3=[2.20,0.67,0.48,0.917];
% P3=[0.0305,0.1171,0.0174,0.0102];
% P3=reshape(P3,[1,1,n3]);
% 
% %The fourth PDF is just the Gaussian matching the first two moments of the
% %original PDF.
% w4=1;
% [mu4,P4]=calcMixtureMoments(mu1,w1,P1);
% 
% DVal=zeros(8,3);
% for curAlg=0:7
%     DVal(curAlg+1,1)=KLDistApproxGaussMix(w1,mu1,P1,w2,mu2,P2,curAlg);
%     DVal(curAlg+1,2)=KLDistApproxGaussMix(w1,mu1,P1,w3,mu3,P3,curAlg);
%     DVal(curAlg+1,3)=KLDistApproxGaussMix(w1,mu1,P1,w4,mu4,P4,curAlg);
% end
% RelErr=abs(bsxfun(@rdivide,bsxfun(@minus,DVal,DVal(1,:)),DVal(1,:)))
% %The rows of RelErr are the algorithms and the PDF chosen for comparison
% %against the original PDF is given by the column. The value obtained
% %through the numeric integration (algorithm 0) is nearly exact, so it is
% %taken as a basis of comparison for the relative error of the
% %approximations. It can be seen that some approximations are much worse
% %than others.
% %We can visualize the above PDFs.
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w1,mu1,P1);
% PDFVals2=GaussianMixtureD.PDF(xVals,w2,mu2,P2);
% PDFVals3=GaussianMixtureD.PDF(xVals,w3,mu3,P3);
% PDFVals4=GaussianMixtureD.PDF(xVals,w4,mu4,P4);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'-r','linewidth',2)
% plot(xVals,PDFVals3,'--g','linewidth',2)
% plot(xVals,PDFVals4,'-.b','linewidth',2)
% legend('10 Elements', '5 Elements', '3 Elements','Matched Gaussian')
%
%EXAMPLE 2:
%Here is an example with bivariare Gaussian mixtures. 
% w1=[0.25;0.5;0.25];
% mu1=zeros(2,2);
% mu1(:,1)=[1;-1];
% mu1(:,2)=[-1;1];
% mu1(:,3)=[0;0];
% P1=zeros(2,2,2);
% P1(:,:,1)=[4/9,  14/45;
%            14/45,4/9];
% P1(:,:,2)=[4/9, 0;
%            0, 4/9];
% P1(:,:,3)=[2/9, -1/9;
%           -1/9, 3/9];
% 
% %The second distribution just throws out the first component.
% w2=w1(2:3);
% w2=w2/sum(w2);
% mu2=mu1(:,2:3);
% P2=P1(:,:,2:3);
% 
% %The third distribution replaces the original distribution with a single
% %moment-matched Gaussian.
% w3=1;
% [mu3,P3]=calcMixtureMoments(mu1,w1,P1);
%       
% DVal=zeros(8,2);
% for curAlg=0:7
%     DVal(curAlg+1,1)=KLDistApproxGaussMix(w1,mu1,P1,w2,mu2,P2,curAlg);
%     DVal(curAlg+1,2)=KLDistApproxGaussMix(w1,mu1,P1,w3,mu3,P3,curAlg);
% end
% RelErr=abs(bsxfun(@rdivide,bsxfun(@minus,DVal,DVal(1,:)),DVal(1,:)))
% %The rows of RelErr are the algorithms and the PDF chosen for comparison
% %against the original PDF is given by the column. The value obtained
% %through the numeric integration (algorithm 0) has a few digits of
% %accuracy, so it is taken as a basis of comparison for the relative
% %error of the approximations. It can be seen that some approximations are
% %much worse than others.
% %We can visualize the above PDFs.
% figure(1)
% clf
% subplot(3,1,1);
% numPoints=250;
% vals=linspace(-3,3,numPoints);
% [X,Y]=meshgrid(vals,vals);
% points=[X(:)';Y(:)'];
% PDFVals=GaussianMixtureD.PDF(points,w1,mu1,P1);
% PDFVals=reshape(PDFVals,numPoints,numPoints);
% surf(X,Y,PDFVals,'EdgeColor','none')
% title('Full PDF')
% %Plot the reduced PDF
% subplot(3,1,2);
% numPoints=250;
% vals=linspace(-3,3,numPoints);
% [X,Y]=meshgrid(vals,vals);
% points=[X(:)';Y(:)'];
% PDFVals=GaussianMixtureD.PDF(points,w2,mu2,P2);
% PDFVals=reshape(PDFVals,numPoints,numPoints);
% surf(X,Y,PDFVals,'EdgeColor','none')
% title('Missing 1 Component')
% %Plot the Gaussian approximation.
% subplot(3,1,3);
% numPoints=250;
% vals=linspace(-3,3,numPoints);
% [X,Y]=meshgrid(vals,vals);
% points=[X(:)';Y(:)'];
% PDFVals=GaussianMixtureD.PDF(points,w3,mu3,P3);
% PDFVals=reshape(PDFVals,numPoints,numPoints);
% surf(X,Y,PDFVals,'EdgeColor','none')
% title('Gaussian Approximation')
%
%REFERENCES:
%[1] J. L. Durrieu, J. P. Thiran, and F. Kelly, "Lower and upper bounds for
%    approximation of Kullback-Leibler divergence between Gaussian mixture
%    models," in Proceedings of the IEEE International Conference on
%    Acoustics, Speech and Signal Processing, Kyoto, Japan, 25-30 Mar.
%    2012, pp. 4833-4836.
%[2] J. Goldberger, S. Gordon, and H. Greenspan, "An efficient image
%    similarity measure based on approximations of KL-divergence between
%    two Gaussian mixtures," in Proceedings of the Ninth IEEE International
%    Conference on Computer Vision, Nice, France, 2003.
%[3] J. R. Hershey and P. A. Olden, "Approximating the Kullback Leibler
%    divergence between Gaussian mixture models," in Proceedings of the
%    IEEE International Conference on Acoustics, Speech and Signal
%    Processing, Honolulu, HI, 15-20 Apr. 2007.
%[4] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(algorithm))
    algorithm=2;
end

d=size(mu1,1);
n1=length(w1);
n2=length(w2);

totalError=[];
exitCode=0;
switch(algorithm)
    case 0%Numeric integration over a bounded region.
        if(nargin>7&&~isempty(param8))
            if(isfield(param8,'numStdDev'))
                numStdDev=param8.numStdDev;
            end
        else
            numStdDev=4.5;
        end

        [mu1Merged,P1Merged]=calcMixtureMoments(mu1,w1,P1);
        [mu2Merged,P2Merged]=calcMixtureMoments(mu2,w2,P2);

        stdDevStep1=numStdDev*sqrt(diag(P1Merged));
        stdDevStep2=numStdDev*sqrt(diag(P2Merged));
        
        P1Inv=zeros(d,d,n1);
        P1InvDet=zeros(n1,1);
        P2Inv=zeros(d,d,n2);
        P2InvDet=zeros(n2,1);
        for k=1:n1
            P1Inv(:,:,k)=inv(P1(:,:,k));
            P1InvDet(k)=det(P1Inv(:,:,k));
        end
        for k=1:n2
            P2Inv(:,:,k)=inv(P2(:,:,k));
            P2InvDet(k)=det(P2Inv(:,:,k));
        end
        
        upperBounds=max(mu1Merged+stdDevStep1,mu2Merged+stdDevStep2);
        lowerBounds=min(mu1Merged-stdDevStep1,mu2Merged-stdDevStep2);

        f=@(x)KLArgGaussMix(x,w1,mu1,P1Inv,P1InvDet,w2,mu2,P2Inv,P2InvDet);

        if(d==1)
            %A different integration routine is used for 1D veruss
            %arbitrary dimensional integrals.
            bounds=[lowerBounds;upperBounds];

            RelTol=[];
            AbsTol=[];
            NPowMin=[];
            NPowMax=[];
            if(nargin>7&&~isempty(param8))
                if(isfield(param8,'RelTol'))
                    RelTol=param8.RelTol;
                end
                
                if(isfield(param8,'AbsTol'))
                    AbsTol=param8.AbsTol;
                end
                
                if(isfield(param8,'NPowMin'))
                    NPowMin=param8.NPowMin;
                end
                                
                if(isfield(param8,'NPowMax'))
                    NPowMax=param8.NPowMax;
                end
            end

            [DVal,totalError,exitCode]=integral1DAdaptiveCC(f,bounds,RelTol,AbsTol,NPowMin,NPowMax);
        else
            maxSearchReg=[];
            AbsTol=[];
            RelTol=[];
            if(nargin>7&&~isempty(param8))
                if(isfield(param8,'RelTol'))
                    RelTol=param8.RelTol;
                end
                
                if(isfield(param8,'AbsTol'))
                    AbsTol=param8.AbsTol;
                end
                
                if(isfield(param8,'maxSearchReg'))
                    maxSearchReg=param8.maxSearchReg;
                end
            end

            [DVal,totalError,exitCode]=integrateUnifCub57Adaptive(f,lowerBounds,upperBounds,maxSearchReg,AbsTol,RelTol);
        end
    case 1%Monte Carlo sampling from Equation 4 of Section 2 of [3].
        if((nargin<8||isempty(param8)))
            numSamp=2e3;
        else
            numSamp=param8;
        end
        
        S1=zeros(d,d,n1);
        S2=zeros(d,d,n2);
        for k=1:n1
            S1(:,:,k)=chol(P1(:,:,k),'lower');
        end
        for k=1:n2
            S2(:,:,k)=chol(P2(:,:,k),'lower');
        end
        
        x1Samp=GaussianMixtureD.randS(numSamp,w1,mu1,S1);
        PDFVals1=GaussianMixtureD.PDFS(x1Samp,w1,mu1,S1);
        
        logVals=log(PDFVals1)-log(GaussianMixtureD.PDFS(x1Samp,w2,mu2,S2));
        %Get rid of points where PDFVals1==0 due to finite precision
        %errors. This probably won't occur. However, eliminating them also
        %gets rid of any unlikely points where both PDFs evaluate to 0 and
        %thus avoids returning NaN values.
        sel=(PDFVals1==0);
        logVals(sel)=0;

        DVal=sum(logVals)/numSamp;
    case 2%The cubature integration approximation from Equation 8 in
          %Section 3 of [3].
        if(nargin<8||isempty(param8))
            [xi,w]=fifthOrderCubPoints(d);
        else
            xi=param8;
        end

        S1=zeros(d,d,n1);
        S2=zeros(d,d,n2);
        for k=1:n1
            S1(:,:,k)=chol(P1(:,:,k),'lower');
        end
        for k=1:n2
            S2(:,:,k)=chol(P2(:,:,k),'lower');
        end
        
        w=w(:).';
        DVal=0;
        for a=1:n1
            %The component over which the expected value is taken.
            xiCur=transformCubPoints(xi,mu1(:,a),S1(:,:,a));

            logVals=log(GaussianMixtureD.PDFS(xiCur,w1,mu1,S1))-log(GaussianMixtureD.PDFS(xiCur,w2,mu2,S2));
            %Get rid of NaN values.
            %logVals(isnan(logVals))=0;
            
            DVal=DVal+w1(a)*sum(logVals.*w);
        end 
    case 3%The variational upper bound from Section 8 of [3].
        psi=zeros(n1,n2);
        phi=zeros(n2,n1);
        
        %We start with psi and phi being the same and equal to what would
        %produce the convexity bound:
        for a=1:n1
            for b=1:n2
                psi(a,b)=w1(a)*w2(b);
                phi(b,a)=psi(a,b);
            end
        end
        
        %Get all pairwise KL-divergences.
        DKLg=zeros(n1,n2);
        for a=1:n1
            for b=1:n2
                DKLg(a,b)=KLDistGauss(mu1(:,a),P1(:,:,a),mu2(:,b),P2(:,:,b));
            end
        end
        expnDKLg=exp(-DKLg);

        numIter=10;
        for curIter=1:numIter
            %Equation 24 for phi.
            for a=1:n1
                denom=sum(psi(a,:).*expnDKLg(a,:));
                for b=1:n2
                    phi(b,a)=w1(a)*psi(a,b)*expnDKLg(a,b)/denom;
                end
            end
        
            %Equation 23 for psi
            for b=1:n2
                denom=sum(phi(b,:));
                for a=1:n1 
                    psi(a,b)=w2(b)*phi(b,a)/denom;
                end
            end
        end
        
        %The first term in Equation 22
        term1=0;
        term2=0;
        for a=1:n1
            for b=1:n2
                val=phi(b,a)*log(phi(b,a)/psi(a,b));
                if(isnan(val))
                    val=0;
                end
                term1=term1+val;

                term2=term2+phi(b,a)*DKLg(a,b);
            end
        end

        DVal=term1+term2;
    case 4%The variational approximation from Section 2.3 of [1] and
          %Section 7 of [3].

        %Compute the necessary KL divergences.
        DKLf=zeros(n1,n1);
        DKLg=zeros(n1,n2);
        for a=1:n1
            DKLf(a,a)=KLDistGauss(mu1(:,a),P1(:,:,a),mu1(:,a),P1(:,:,a));
            for ap=(a+1):n1
                DKLf(a,ap)=KLDistGauss(mu1(:,a),P1(:,:,a),mu1(:,ap),P1(:,:,ap));
                DKLf(ap,a)=KLDistGauss(mu1(:,ap),P1(:,:,ap),mu1(:,a),P1(:,:,a));
            end

            for b=1:n2
                DKLg(a,b)=KLDistGauss(mu1(:,a),P1(:,:,a),mu2(:,b),P2(:,:,b));
            end
        end

        %Evaluate the sum in Equation 18 of [1], which is the same as
        %Equation 20 of [3].
        DVal=0;
        w1=w1(:).';
        w2=w2(:).';
        for a=1:n1
            DVal=DVal+w1(a)*(log(sum(w1.*exp(-DKLf(a,:))))-log(sum(w2.*exp(-DKLg(a,:)))));
        end
    case 5%The product of Gaussians approximation in [1] and Section 5 of
          %[3].
          
        P1Inv=zeros(d,d,n1);
        P2Inv=zeros(d,d,n2);
        for k=1:n1
            P1Inv(:,:,k)=inv(P1(:,:,k));
        end
        for k=1:n2
            P2Inv(:,:,k)=inv(P2(:,:,k));
        end

        %We have to first compute the z and t terms. The logarithms of the z
        %and t terms are given in Equation 22 in [1]. However, that
        %expression is incorrect; the corrected expression is used here

        zaap=zeros(n1,n1);
        tab=zeros(n1,n2);
        for a=1:n1
            for ap=a:n1
                diff=mu1(:,a)-mu1(:,ap);
                Sigmad=inv(P1Inv(:,:,a)+P1Inv(:,:,ap));
                const=sqrt(det(2*pi*Sigmad))/(sqrt(det(2*pi*P1(:,:,a))*det(2*pi*P1(:,:,ap))));
                
                B=P1Inv(:,:,a)*Sigmad*P1Inv(:,:,ap);
                zaap(a,ap)=const*exp(-(1/2)*diff'*B*diff);
                zaap(ap,a)=zaap(a,ap);
            end

            for b=1:n2
                diff=mu1(:,a)-mu2(:,b);
                Sigmad=inv(P1Inv(:,:,a)+P2Inv(:,:,b));
                const=sqrt(det(2*pi*Sigmad))/(sqrt(det(2*pi*P1(:,:,a))*det(2*pi*P2(:,:,b))));
                
                B=P1Inv(:,:,a)*Sigmad*P2Inv(:,:,b);
                tab(a,b)=const*exp(-(1/2)*diff'*B*diff);
            end
        end
        
        %Evaluate the sum in Equation 12 of [1], which is the same as
        %Equation 12 of [3].
        DVal=0;
        w1=w1(:).';
        w2=w2(:).';
        for a=1:n1
            DVal=DVal+w1(a)*(log(sum(w1.*zaap(a,:)))-log(sum(w2.*tab(a,:))));
        end
    case 6
        DProd=KLDistApproxGaussMix(w1,mu1,P1,w2,mu2,P2,5);
        DVar=KLDistApproxGaussMix(w1,mu1,P1,w2,mu2,P2,4);
        DVal=(DProd+DVar)/2;
    case 7%Goldberger's approximation from Equation 2 of [2].
        DKLg=zeros(n1,n2);
        for a=1:n1
            for b=1:n2
                DKLg(a,b)=KLDistGauss(mu1(:,a),P1(:,:,a),mu2(:,b),P2(:,:,b));
            end
        end
        
        w2=w2(:).';

        %Perform the minimization in Equation 1 of [2], which is also given
        %in Equation 14 of [3] and simultaneously evaluate Equation 2 of
        %[2], which is the same as Equation 15 in [3].
        DVal=0;
        logw2=log(w2);
        for a=1:n1
            minVal=min(DKLg(a,:)-logw2);
            DVal=DVal+w1(a)*(log(w1(a))+minVal);
        end
    otherwise
        error('Unknown algorithm specified.')
end
end

function val=KLArgGaussMix(x,w1,mu1,P1Inv,P1InvDet,w2,mu2,P2Inv,P2InvDet)
%%KLSAMPGAUSSMIX The KL divergence is an integral. This evaluates the
%                argument of the integral for the KL divergences between
%                two Gaussian mixtures.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

PDF1=GaussianMixtureD.PDFI(x,w1,mu1,P1Inv,P1InvDet);
PDF2=GaussianMixtureD.PDFI(x,w2,mu2,P2Inv,P2InvDet);

val=PDF1.*(log(PDF1)-log(PDF2));
val(isnan(val))=0;

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
