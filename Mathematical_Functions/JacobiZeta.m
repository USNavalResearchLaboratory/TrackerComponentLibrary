function intVal=JacobiZeta(beta,m)
%%JACOBIZETA Evaluate Jacobi's zeta function. This is defined to be
%            Z(beta,m)=E(beta|m)-E(pi/2,m)/F(pi/2,m)*F(beta|m), where F is
%            the incomplete elliptic integral fo the first kind and E is
%            the incomplete elliptic integral of the second kind. The
%            Jacobi zeta function is often denoted b Z. The function arises
%            in various obscure scenarios, such as in [2].
%
%INPUTS: beta, m The two real or complex parameters of the Jacobi zeta
%                function.
%
%OUTPUTS: intVal The scalar value fo the jacobi zeta function with the
%                given parameters.
%
%This function implements the formula given in Equation 63 of [1].
%
%REFERENCES:
%[1] B. C. Carlson, "Numerical computation of real or complex elliptic
%    integrals," Numerical Algorithms, vol. 10, no. 1, pp. 13-26, 1995.
%[2] A. J. Brizard, "Jacobi zeta function and action-angle coordinates
%    for the pendulum," Communications in Nonlinear Science and
%    Numerical Simulation, vol. 18, no. 3, pp. 511-518, Mar. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    sinBeta=sin(beta);
    K=ellipIntInc1Kind(pi/2,m);
    intVal=(m/3)*sinBeta*cos(beta)*sqrt(1-m*sinBeta^2)*symIntThirdKind(0,1-m,1,1-m*sinBeta^2)/K;
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
