function val=BesselI1RatioApprox(kappa,algorithm)
%%BESSELI1RATIOAPPROX Approximate the ratio of modified Bessel functions of
%              the first kind of the form x=I_{1}(kappa)/I_{0}(kappa).
%              For a more exact solution, use the function BesseliRatio,
%              which can also handle different subscripts (orders) for the
%              modified Bessel functions of the first kind.
%
%INPUTS: kappa The real argument of the modified Bessel function of the
%              first kind kappa>=0.
%    algorithm A parameter selecting the approximation that is used.
%              Possible values are:
%              0 (The default if omitted or an empty matrix is passed) Use
%                the approximation used in Equations 8 and 10 of [1].
%                Equation 8 in [1] is the first two terms of an expansion
%                given on page 290 of [2].
%              1 Use the approximations in Equations 5 and 7 of [3].
%                Equation 5 of [3] is related to a Taylor series expansion
%                given on page 289 of [2].
%
%OUTPUTS: x An approximation to the ratio I_{1}(kappa)/I_{0}(kappa).
%
%REFERENCES:
%[1] G. Stienne, S. Reboul, M. Azmani, J. B. Choquel, and M. Benjelloun, A
%    multi-temporal multi-sensor circular fusion filter," Information
%    Fusion, vol. 18, pp. 86-100, Jul. 2014.
%[2] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%[3] G. Stienne, S. Reboul, J. B. Choquel, and M. Benjelloun, "Circular
%    data processing tools applied to a phase open loop architecture for
%    multi-channel signals tracking," in Position Location and Navigation
%    Symposium, Myrtle Beach, SC, 23-26 Apr. 2012, pp. 633-642.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(algorithm))
        algorithm=0; 
    end
    
    switch(algorithm)
        case 0%The method of [1].
            if(kappa>=0.6)
                val=1-1/(2*kappa);
            else
                val=exp(-1/(2*kappa));
            end
        case 1
            if(kappa>=0.6)
                val=(1-3/(8*kappa)-15/(128*kappa^2))/(1+1/(8*kappa)+9/(128*kappa^2));
            else
                 val=kappa/2;
            end
        otherwise
            error('Unknown Algorithm Specified.')
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
