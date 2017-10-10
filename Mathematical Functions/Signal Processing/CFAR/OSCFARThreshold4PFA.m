function T=OSCFARThreshold4PFA(PFA,N,k)
%%OSCFARTHRESHOLD4PFA Determine the threshold needed to obtain a chosen
%           probability of false alarm in uniform clutter when using the
%           order-statistics constant false alarm (OS-CFAR) detector with a
%           simple exponential noise model as used in [1].
%
%INPUTS: PFA The desired probability of false alarm (0<PFA<1) to use for
%          the CFAR detector. A matrix of values can be passed in which
%          case the output is a matrix evaluated at each point.
%        N The total number of cells that contribute to the test region.
%          These are usually cells that are around the point being tested,
%          some distance outside of a guard region.
%        k The integer order to use k. That is, the kth largest sample in
%          the test region is used as the test statistic.
%
%OUTPUTS: T The threshold for the given parameters, T>0.
%
%This function is based on Equation 37 of [1] with S=0. We subtract the
%right-hand side so that we have PFA-prod...=0. Then, we multiple by the
%inverse of the product so that we have PFA*prod...-1=0. The resulting
%expression is a polynomial. We just find the roots of the polynomial and
%then throw out all solutions that are obviously wrong (complex, negative,
%etc). Then if, there are any solutions left, we take the one that produces
%a PFA closes to the desired value.
%
%The function OSCFARPFA4Threshold is the inverse of this function.
%
%EXAMPLE:
%This example is notable, because with N and k so large, many other
%techniques would fail to produce a result.
% N=180;
% k=179;
% PFA=1e-8;
% T=OSCFARThreshold4PFA(PFA,N,k)
%
%REFERENCES
%[1] P. P. Gandhi and S. A. Kassam, "Analysis of CFAR processors in
%    nonhomogeneous background," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 24, no. 4, pp. 427-445, Jul. 1988.
%
%February 2017  David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVals=numel(PFA);
T=zeros(size(PFA));

for curVal=1:numVals
    PFACur=PFA(curVal);

    %Determine the product of the numerators
    prodPoly=1;
    for i=0:(k-1)
        %Polynomial multiplication
        prodPoly=conv(prodPoly,[1/(N-i),1]);
    end
    prodPoly=PFACur*prodPoly;
    prodPoly(end)=prodPoly(end)-1;

    %Solve the polynomial.
    rootVals=roots(prodPoly);
    %Throw out all negative and zeor roots.
    rootVals(rootVals>0);

    %Throw out all roots that are obviously imaginary and only keep the real
    %part of the remaining ones.
    sel=imag(rootVals)<eps(real(rootVals));
    rootVals=real(rootVals(sel));

    numRoots=numel(rootVals);
    if(numRoots>1)
        %Choose the one that is closest to the desired value.
        [~,idx]=min(abs(OSCFARPFA4Threshold(rootVals,N,k)-PFACur));
        T(curVal)=rootVals(idx);
    elseif(numRoots==0)
        error('Finite precision errors prevent a solution from being found')
    else
        T(curVal)=rootVals;
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
