function [J1Idx,J2Idx]=sampleStratonovichIntegrals(m,deltaT,p)
%%SAMPLESTRATONOVICHINTEGRALS Return a random sample of Stratonovich
%         integrals with 1 and 2 subscripts. These are multiple
%         Stratonovich integrals of 1 (not of a general function). The
%         assumed Wiener process can be multidimensional, in which case
%         this function returns the necessary cross terms too. All values
%         are generated with the same assumed Wiener process. Compare this
%         to sampleItoIntegrals.
%
%INPUTS: m The scalar dimensionality of the Wiener process, m>=1.
%   deltaT The time interval over which the integrals are taken.
%        p The number of terms to use in the approximation. The default if
%          omitted or an empty matrix is passed is 100.
%
%OUTPUTS: J1Idx A structure whose members are all Stratonovich integrals
%               with 1 subscript. Members are:
%               J0 for J_0 (scalar) This is the deterministic integral of
%                                   the region 0->deltaT.
%               Jj for J_j (mX1)
%         J2Idx A structure whose members are all Ito integrals with 2
%               subscripts. Members are:
%               J00 for J_{0,0} (scalar) This is the deterministic second
%                                        integral of the region 0->deltaT.
%               J0j for J_{0,j} (mX1)
%               Jj0 for J_{j,0} (mX1)
%               Jj1j2 for J_{j1,j2} (mXm)
%
%This implements the approximations of Chapter 5.8 of [1]. 
%
%The Stratonovich integral is defined in a manner similar to a
%Riemann-Stieljes integral. That is, the integral of some function f(x),
%when taking N discrete steps, is
%sum_{j=0}^{N-1} f(tau_j)*(W_{t_{j+1}}-W_{t_j})
%where W_{t_j} is the value of a Wiener process at discrete time t_j and
%tau_j=(t_{j+1}+t_j)/2. If one wished to evaluate these integrals,
%brute-force, one would have to simulate 2*p-1 values of the Wiener process
%in order to be able to get values between the steps in the double
%integrals.
%
%EXAMPLE:
%Here, we demonstrate how some of the second moments, when converted to the
%equivalent Ito formulation, match the value given by ItoIntegral2ndMoment.
%The formula of Equation 2.36 in Chapter 5.2 of [1] is used to convert
%Jj1j2 into the equivalent Itô form.
% m=2;
% deltaT=2;
% p=50;
% numMonteCarlo=1e5;
% Ij1j2Moment=zeros(m,m,numMonteCarlo);
% for curRun=1:numMonteCarlo
%     [J1Idx,J2Idx]=sampleStratonovichIntegrals(m,deltaT,p);
%     J0=J1Idx.J0;
% 
%     Jj1j2=J2Idx.Jj1j2;
%     %Equation 2.36 in Chapter 5.2
%     Ij1j2=Jj1j2-(1/2)*eye(m,m)*J0;
% 
%     Ij1j2Moment(:,:,curRun)=Ij1j2.^2;
% end
% meanSim=mean(Ij1j2Moment,3)
% exactSol=ItoIntegral2ndMoment([1,1],[1,1],deltaT)
%The exact solution is 2 (It is the same for all elements) and the meanSim
%will have values that are usually within 0.03 of 2.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(p))
    p=100;
end

sqrtDeltaT=sqrt(deltaT);
xi=randn(m,1);

J0=deltaT;
Jj=sqrtDeltaT*xi;

J1Idx.J0=J0;
J1Idx.Jj=Jj;

if(nargout>1)
    zeta=randn(m,p);
    eta=randn(m,p);
    mup=randn(m,1);
    r=1:p;
    rhop=1/12-(1/(2*pi^2))*sum(1./r.^2);

    J00=deltaT^2/2;

    a0=-1/pi*sqrt(2*deltaT)*sum(bsxfun(@times,1./r,zeta),2)-2*sqrt(deltaT*rhop)*mup;
    Jj0=(deltaT/2)*(sqrtDeltaT*xi+a0);
    J0j=(deltaT/2)*(sqrtDeltaT*xi-a0);

    Jj1j2=zeros(m,m);
    A=zeros(m,m);
    for j1=1:m
        for j2=1:m

            A(j1,j2)=0;
            for k=1:p
                A(j1,j2)=A(j1,j2)+(1/k)*(zeta(j1,k)*eta(j2,k)-eta(j1,k)*zeta(j2,k)); 
            end
            A(j1,j2)=A(j1,j2)/(2*pi);

            Jj1j2(j1,j2)=(1/2)*deltaT*xi(j1)*xi(j2)-sqrtDeltaT/2*(a0(j2)*xi(j1)-a0(j1)*xi(j2))+deltaT*A(j1,j2);
        end
    end

    J2Idx.J00=J00;
    J2Idx.Jj0=Jj0;
    J2Idx.J0j=J0j;
    J2Idx.Jj1j2=Jj1j2;
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
