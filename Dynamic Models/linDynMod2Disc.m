function [F,Q,u]=linDynMod2Disc(T,A,D,alpha)
%%LINDYNMOD2DISC Given a continuous-time linear dynamic model with constant
%                drift A and diffusion D matrices and a constant input
%                alpha such that dx=(A*x+u)dt+D*dW, where x is the state,
%                (so dx is the differential state) W is a Wiener process,
%                find the equivalent state transition matrix F, process noise
%                covariance matrix and discrete time input u for a
%                prediction over an interval of duration T. The
%                discrete-time propagation equation is x_{t+T}=F*x_t+u+v,
%                where v is zero-mean Gaussian noise with covariance matrix
%                Q.
%
%INPUTS: T The duration of the prediction interval over which the
%          continuous-time system should be discretized.
%        A The xDimXxDim constant drift matrix, which cannot be singular if
%          a control input u is provided.
%        D The xDimXwDim constant diffusion matrix.
%    alpha The constant continuous-time xDimX1 (control) input. If omitted,
%          it is assumed there is no control input.
%
%OUTPUTS: F The xDimXxDim state transition matrix of the discretized
%           system.
%         Q The xDimXxDim process noise covariance matrix of the
%           discretized system.
%         u The xDimX1 discrete-time control input.
%
%Chapter 6.2.2 of [1] describes how such discretization is performed. The F
%matrix is easily found using a matrix exponential. However, the Q matrix
%requires an integral over a matrix exponential. In Problem 9.3.4 of
%Chapter 9 of [2], a technique for evaluating such a matrix exponential is
%presented. That method is used here. It originally comes from [3]. The
%control input solution is standard from (non-stochastic) linear systems
%theory as given in Chapter 4.2 of [4]. The non-stochastic solution can be
%applied to the stochastic system, because in a linear system, the effects
%of the noise can just be superpositioned on the effects of the non-
%stochastic system.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%[3] C. F. Van Loan, "Computing Integrals Involving the Matrix
%    Exponential," IEEE Transactions on Automatic Control, vol. AC-23, no.
%    3, pp. 395-404, Jun. 1967.
%[4] C.-T. Chen, Linear System Theory and Design, 4th ed. New York: Oxford
%    University Press, 2013.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(A,1);

%The state transition matrix is straightforward and is in many texts. It is
%the same as in the noise-free case.
F=expm(A*T);

%The integral for the process noise covariance matrix is evaluated using
%the track in the Golub and Van Loan book.
P=D*D';
M=[-A,              P;
   zeros(xDim,xDim),A'];
FExp=expm(M*T);%Matrix exponential.

F12=FExp(1:xDim,(xDim+1):end);
F22=FExp((xDim+1):end,(xDim+1):end);
Q=F22'*F12;

%The control input is just the integral from 0 to T of exp(A*t)*alpha dt.
if(nargin>3)
    u=matExpInt(T,A,alpha);
else%If no control input is given.
    u=zeros(xDim,1);
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
