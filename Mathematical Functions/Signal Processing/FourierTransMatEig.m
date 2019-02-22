function [E,lambda,Fa]=FourierTransMatEig(N,order,a)
%%FOURIERTRANSMATEIG Obtain the eigenvectors of a matrix to perform the
%                 discrete Fourier transform. Optionally one can obtain the
%                 eigenvalues of the matrix anf the actual matrix raised to
%                 a desired power. The discrete Fourier transform as
%                 obtained using fft(x) is given by Fa*x*sqrt(N) if x is
%                 length N and one obtains Fa via
%                 [~,~,Fa]=FourierTransMatEig(N,[],a);
%                 Though one could obtain a Fourier matrix with the command
%                 Fa=fft(eye(N,N))/sqrt(N), one would quickly find that the
%                 eigenvalue decomposition is not as numerically reliable
%                 as desired and the accuracy of matrix powers is also not
%                 as high as desired. Hence the reason for this function.
%
%INPUTS: N The length of the sequence whose Fourier transform is desired.
%          If the matrix is too large, overflow errors will occur.
%    order The matrix is approximated using an expansion to the selected
%          order. If this parameter is omitted, the default value of N/2 is
%          used. Order can range from 2 to N-1.
%        a If desired, the discrete Fourier transform matrix can be raised
%          to a power. This specified the power. If omitted or an empty
%          matrix is passed, the default of 1 is used.
%
%OUTPUTS: E The NXN matrix of eigenvectors of the discrete Fourier
%           transform matrix.
%    lambda The NX1 set of eigenvalues of the discrete Fourier transform
%           matrix.
%        Fa The discrete Fourier transform matrix raised to the power a.
%
%This implements the algorithm described in Section 6 of [1].
%
%For a=1, Fa is more accurately determined using Fa=fft(eye(N,N))/sqrt(N).
%However, the eigenvectors and eigenvalues obtained using
%[E,Lambda]=eig(Fa) will be less accurate. 
%
%REFERENCES:
%[1] A. Bultheel and H. E. Martínez Sulbaran, "Computation of the
%    fractional Fourier transform," Applied and Computational Harmonic
%    Analysis, vol. 16, no. 3, pp. 182-202, May 2004.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(a))
    a=1;
end

if(nargin<2||isempty(order))
    order=N/2;
end
order=min(max(2,order),N-1);

%%STEP 1: Construct matrix H from Equation 6

%This vector is repeatedly used to compute the coefficients of d_p(x). It
%is defined on page 193. The use of d_p to get D_p is mentioned before
%Equation 6.
dVec=[1;-2;1];
dp=1;%d_0(x)
polySum = 0;
polyTemp=zeros(1,N);

%Rather than directly compute D_p and C_p, we compute the sum of the
%polynomials that make up D_p in equation six. In other words, we are
%essentially computing the sum of all of the D_p plus
coeff=1/2;
for p=1:(order/2)
    %Update the coefficient that is outside of D_p+C_p in Equation 6.
    if(p>1)
        coeff=coeff*(1-2*p+p^2)/(2*p-4*p^2);
    end

    %The current polynomial term to add to the D_p sum.
    dp=conv(dVec,dp);
    
    %The indexation is consistent with the ordering used
    polyTemp([N-p+1:N,1:(p+1)])=dp;
    polySum=polySum+coeff*polyTemp;
end
%As noted, the c_0 term, which is the constant term in d_p is removed.
%Here, we remove it for the whole sum.
polySum(1)=0;

polyTemp=[polySum(N:-1:2).';polySum(:)];
idx=bsxfun(@plus,(0:(N-1))',repmat(N:-1:1,N,1));
Cp=polyTemp(idx);
Dp=diag(real(fft(polySum)));
H=Cp+Dp;

%%STEP 2: Construct the transformation matrix V from Equations 3 and 4

if(mod(N,2)==1)%If odd
    %Use Equation 3
    r=(N-1)/2;
    Ir=eye(r,r);
    Jr=fliplr(Ir);
    
    V=(1/sqrt(2))*[sqrt(2),zeros(1,2*r);
                   zeros(2*r,1),[Ir,Jr;
                                 Jr,-Ir]];
else
    %Use Equation 4. Note that the 1 should have been a sqrt(2).
    r=(N-2)/2; 
    Ir=eye(r,r);
    Jr=fliplr(Ir);
    V=(1/sqrt(2))*[sqrt(2),zeros(1,2*r+1);
                   zeros(2*r+1,1),[Ir,        zeros(r,1),   Jr;
                                   zeros(1,r),sqrt(2),            zeros(1,r);
                                   Jr,        zeros(r,1),   -Ir]];
end

%%STEP 3: Compute Blocks Ev and Od, defined before Equation 3
r=fix(N/2);
VHV=V*H*V';
Ev=VHV(1:r+1,1:r+1);
Od=VHV((r+2):N,(r+2):N);

%STEP 4: Find the eigenvectors of Ev and the Eigenvectors of Od. THey must
%be sorted.
[Ve,eigE]=eig(Ev);
[~,idx] = sort(diag(eigE),'descend');
Ve=Ve(:,idx);

[Vo,eigD]=eig(Od); 
[~,idx]=sort(diag(eigD),'descend');
Vo=Vo(:,idx);

%STEP 5: Transform the eigenvectors to get E as in Equation 5.
E=zeros(N,N);
E(1:r+1,1:r+1)=Ve;
E(r+2:N,r+2:N)=Vo;
E=V*E;

%STEP 6: Interlace the eigenvectors.

idx=[1:r+1;r+2:2*r+2];
idx=idx(:);
if(mod(N,2)==1)
    idx(N+1)=[];
else
    idx([N,N+2])=[];
end

E=E(:,idx);

%STEP 7
if(nargout>1)
    isEven=(mod(N,2)==0);
    lambda=exp(-1j*pi/2*a*([0:(N-2),N-1+isEven])).';
end

if(nargout>2)
    %The following is the same as Fa=E*diag(Lambda)*E';
    Fa=E*(lambda.*E');
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
