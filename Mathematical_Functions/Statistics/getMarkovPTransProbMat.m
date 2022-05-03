function Lambda=getMarkovPTransProbMat(A,T)
%%GETMARKOVPTRANSPROBMAT Get the transition probability matrix for a
%           continuous-time Markov process given given the transition
%           density matrix for the process and the time delay. This can be
%           useful in certain implementations of the interacting multiple
%           model (IMM) filter that have variable revisit rates.
%
%INPUTS: A The transition density of the continuous-time Markov process.
%          The sum of the elements on every row must be zero and all of
%          the diagonal elements must be negative. The values of the
%          diagonal elements can be linked to the mean sojourn time in a
%          given state, as discussed below.
%        T The time interval over which the transition probability matrix
%          is taken to consider.
%
%OUTPUTS: Lambda The transition probability matrix over the given time
%                interval. Lambda(i,j) is the probability that if the
%                system was originally in state i, that after a time
%                duration of T it will be in state j.
%
%As discussed in Chapter 16.2 of [1], the negative inverses of the diagonal
%elements of A are the mean sojourn times in each state. Given that one is
%in state i at time zero, the probability of leaving state i during a time-
%duration T is 1-exp(-A(i,i)*T) as per Equation 16-10 of Papoulis. The off-
%diagonal terms in A affect to which state a transition will be made when
%switching.
%
%The transformation from the matrix A given T to a transition probability
%matrix comes from the solution to the Kolmogorov equations given in
%Chapter 16.2 of the above reference. Specifically, from Equation 16-29.
%
%EXAMPLE:
%Consider two states. 
% A=[-0.05, 0.05;
%     0.1, -0.1];
% T=1;
% Lambda=getMarkovPTransProbMat(A,T);
%The resulting matrix is thus
%Lambda=[0.9536    0.0464;
%         0.0929    0.9071];
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Equation 16-27.
    Lambda=expm(A*T);
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
