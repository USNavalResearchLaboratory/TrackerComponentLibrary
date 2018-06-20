function S=impedMat2ScatterMat(Z,zPortImped)
%%IMPEDMAT2SCATTERMAT Convert an impedance matrix to a scattering matrix
%           for an N-port network. Impedance matrices relate voltage and
%           current as V=Z*I, where V is an NX1 vector representing the
%           time-harmonic voltages at n ports in a network and I is an NX1
%           vector representing the currents at N ports in a network. The
%           scattering matrix is of the form  b=S*a, where a represents
%           "power waves" going into the network and b represents power
%           waves coming out of the network. Power waves, having units of
%           sqrt(Watts), are discussed in more detail below.
%
%INPUTS: Z The NXN impedance matrix of the network. The units will
%           typically be Ohms.
% zPortImped The NX1 or 1XN vector of the port impedances that define the
%          scattering matrix. The units are typically Ohms. These can be
%          complex. If this parameter is omitted or an empty matrix is
%          passed, then the default of 50 Ohms is used. If a scalar is
%          passed, then all ports are set to that impedance.
%
%OUTPUTS: S The NXN complex scattering matrix corresponding to Z. The
%          entries of this are unitless. The value in S(i,j) is b(i)/a(j),
%          where b is a length N vector of power wave amplitudes exiting
%          the ports in the network and a is a length N vector of the power
%          wave amplitudes entering the ports. The amplitudes, with units
%          sqrt(Watts) are defined below.
%
%Power waves and the scattering matrix are described in [1]. A
%disambiguation between power waves and other waves, such as pseudowaves,
%is given in [2]. The definition of power waves used here is as in [1]:
% a(i)=(V(i)+zPortImped(i)*I(i))/(2*sqrt(abs(real(zPortImped(i)))))
% b(i)=(V(i)-conj(zPortImped(i))*I(i))/(2*sqrt(abs(real(zPortImped(i)))))
%where V(i) and I(i) are the Voltage and current at the ith port. If
%zPortImped is the actual set of impedances to which the ports of the
%network are attached, then the power going into the ith port is
%P=sign(real(zPortImped(i)))*(abs(a(i))-abs(b(i)))
%The inverse transformation from V and I to power waves a and b is
%V(i)=(sign(real(zPortImped(i)))/sqrt(abs(real(zPortImped(i)))))*(conj(zPortImped(i))*a(i)+zPortImped(i)*b(i))
%I(i)=(sign(real(zPortImped(i)))/sqrt(abs(real(zPortImped(i)))))*(a(i)-b(i))
%Note that some books use other definitions and often times zPortImped is
%assumed to be all real, often a fixed 50 Ohms.
%
%The ith element of the impedance matrix is V(i)/I(i). The ith element of
%the scattering matrix is b(i)/a(i). The relation between the impedance
%matrix used here is derived in a manner similar to that used in [3], which
%uses a less general definition for the port impedances.
%
%REFERENCES:
%[1] K. Kurokawa, "Power waves and the scattering matrix," IEEE
%    Transactions on Microwave Theory and Techniques, vol. 13, no. 2, pp.
%    194-202, Mar. 1965.
%[2] R. B. Marks and D. F. Williams, "A general waveguide circuit theory,"
%    Journal of the National Institute of Standards and Technology, vol.
%    97, no. 5, pp. 533-562, Sep.-Oct. 1992.
%[3] R. A. Speciale, "Derivation of the generalized scattering parameter
%    renormalization transformation," in Proceedings of the IEEE
%    International Symposium on Circuits and Systems, vol. 1, Houston, TX,
%    28-30 Apr. 1980, pp. 166-169.
%
%September 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(zPortImped))
    zPortImped=50; 
end

N=size(Z,1);

if(isscalar(zPortImped))
   zPortImped=zPortImped*ones(N,1); 
end

sqrtRezPort=sqrt(abs(real(zPortImped)));
Z0sqrtR=diag(sqrtRezPort);
Z0sqrtRInv=diag(1./sqrtRezPort);

S=Z0sqrtRInv*((Z-diag(conj(zPortImped)))/(Z+diag(zPortImped)))*Z0sqrtR;

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
