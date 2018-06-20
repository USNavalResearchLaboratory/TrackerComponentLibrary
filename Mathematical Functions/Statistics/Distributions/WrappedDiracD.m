classdef WrappedDiracD
%%WRAPPEDDIRACD Functions to handle the wrapped Dirac distribution, which
%           is a circular distritbuion (unique over a value range of 2*pi).
%Implemented methods are: CDF, trigMoment, params4TrigMoment, rand
%
%The wrapped Dirac distribution plays a role in some filtering algorithms
%on circular domains, such as in [1].
%
%REFERENCES:
%[1] G. Kurz, I. Gilitschenski, and U. D. Hanebeck, "Recursive nonlinear
%    filtering for angular data based on circular distributions," in
%    American Control Conference, Washington, DC, 17-19 Jun. 2013.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=CDF(x,d,w,startingPoint)
%%CDF Evaluate the cumulative distribution function of a wrapped Dirac
%     distribution having the given parameters.
% 
%INPUTS: x A vector or matrix of points at which one wishes to evaluate the
%          wrapped normal distribution.  
%        d The 1XnumComp or numCompX1 vector of  real wrapped Dirac points.
%        w The numCompX1 weights associated with the wrapped Dirac points
%          sum(w)=1; w>=0.
%
%OUTPUTS: vals The values of the wrapped normal CDF evaluated at the points
%              in x. This has the same dimensionality as x.
%
%This function  basically just sums weights of the points that are less
%than x when starting at startingPoint. Note that the CDF of a circular
%distirbution is periodic.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
   
    if(nargin<4||isempty(startingPoint))
       startingPoint=0; 
    end
    
    %Shift so that startingPoint is at the beginning of the inverval.
    d=mod(d,2*pi);%Get d in the range of 0 to 2*pi.
    
    d=d+(d<startingPoint)*2*pi;
    x=mod(x,2*pi);
    %Shift the points.
    x=x+(x<startingPoint)*2*pi;
    
    %The points in d do not need to be sorted in any particular manner.
    val=arrayfun(@(y)sum(w(d<=y)),x);
end

function [rho,theta,R]=trigMoment(n,d,w)
%%TRIGMOMENT Compute the nth trigonometric moment of a wrapped Dirac
%            distribution having the given parameters. A direction theta on
%            the unit circle can be modeled as a complex quantity having
%            unit magnitude as exp(1j*theta). The nth raw complex
%            trigonometric moment is the expected value of exp(1j*n*theta),
%            where the integral for the expected value is taken over any
%            interval of length 2*pi.
%
%INPUTS: n The order of the moment desired. This is >=1.
%        d The 1XnumComp or numCompX1 vector of  real wrapped Dirac points.
%        w The numCompX1 weights associated with the wrapped Dirac points
%          sum(w)=1; w>=0.
%
%OUTPUTS: rho The (complex) mean resultant value for the nth moment. Note
%             that rho=R*exp(1j*theta).
%       theta The (real) trigonometric mean angle in radians for the nth
%             moment. This is between -pi and pi.
%           R The (real) mean resultant length for the nth moment.
%
%This method just calls findTrigMomentFromSamp.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    rho=sum(exp(1i*n*d(:)).*w(:));
    if(nargout>1)
        theta=angle(rho);
        R=abs(rho);
    end
end

function [d,w]=params4TrigMoment(rhoList,numComp,lambda)
%%APPROXCIRCDISTBYDIRAC Approximate a circular distribution by a
%           wrapped Dirac mixture with 2, 3, or 5 components such that the
%           wrapped Dirac mixture matches the first moment and in the case
%           of 5 components, the magnitude but not the phase of the second
%           moment. If the distribution is symmetric, then the phase of the
%           fit with both momenets should be exact.
%
%INPUTS: rhoList A length-1 (or length-2) vector of the first (and second)
%                complex trigonometric moment of the circular distribution
%                that is to be approximated. These are the complex mean
%                resultant lengths. The absolute value of each entry in
%                rhoList(1) must be <=1. See below for additional
%                limitation on rhoList(2).
%        numComp The number of Dirac components desired in the
%                approximation. This can be 2, 3, or 5. If 5, then rhoList
%                must be length 2. If less than 5, then rhoList must be
%                length 1. The default if this parameter is omitted or an
%                empty matrix is passed is 2 if rhoList is length 1 and 5
%                if rhoList is length 2. Note that the PHASE of rhoList(2)
%                is not matched.
%         lambda This parameter is only used if numComp==5. This is a value
%                between 0 and 1 that affects how much weight is given to a
%                point at the mean. The default if this parameter is
%                omitted is 0.8. Note that this value is not equal to the
%                weight on the point at the origin.
%
%OUTPUTS: d The 1XnumComp locations of the wrapped Dirac points. These
%           range from -pi to pi.
%         w The numCompX1 weights associated with the wrapped Dirac points
%           sum(w)=1; w>=0.
%
%Approximations for 2 and 3 components are given in [1] and [2]. The
%approximation for 5 components is Algorithm 1 in [2].
%
%In the case of five components, when trying to match rhoList(2), a
%solution does not always exist. As noted in [3], a solution only exists if
%the matrix 
% X=[1,           rhoList(1),     rhoList(2);
%  rhoList(1)',   1,              rhoList(1);
%  rhoList(2)',   rhoList(1)',    1];
%is postivie semidefinite. This is not noted in [1] and [2]. This function
%tests whether that matrix is positive semidefinite and raises and error if
%it is not. However, rather than actually using rhoList(2), the test is
%done on the matrix with the phase of rhoList(2) optimized to provide the
%maximum determininant. This is because this function only matches the
%magnitude, not the phase of rhoList(2).
%
%EXAMPLE:
%Here, we demonstrate that we can match moments using the params4TrigMoment
%and get back the same values using the trigMoment function.
% rhoList(1)=0.3+0.2*1j;
% rhoList(2)=0.1;
% %We fit two points to match the first moment
% [d2,w2]=WrappedDiracD.params4TrigMoment(rhoList,2);
% rho2pt=WrappedDiracD.trigMoment(1,d2,w2)
% %and find that rho2pt is rhoList(1)
% 
% %We fit three points to match the first moment.
% [d3,w3]=WrappedDiracD.params4TrigMoment(rhoList,3);
% rho3pt=WrappedDiracD.trigMoment(1,d3,w3)
% %and find that rho3pt is rhoList(1)
% 
% %We fit five points to match the first moment and the magnitude of the
% %second moment.
% [d5,w5]=WrappedDiracD.params4TrigMoment(rhoList,5,1);
% rho5pt=WrappedDiracD.trigMoment(1,d5,w5)
% rho5pt2=WrappedDiracD.trigMoment(2,d5,w5)
% abs(rho5pt2)
% %We find that rho5pt matches rhoList(1) and rho5pt2 has the same magnitude
% %as rhoList(2), but a different complex phase.
%   
%REFERENCES:
%[1] G. Kurz, I. Gilitschenski, and U. D. Hanebeck, "Deterministic
%    approximation of circular densities with symmetric Dirac mixtures
%    based on two circular moments," in The 17th International Conference
%    on Information Fusion, Salamanca, Spain, 7-10 Jul. 2014.
%[2] G. Kurz, I. Gilitschenski, R. Y. Siegwart, and U. D. Hanebeck,
%    "Methods for deterministic approximation of circular densities,"
%    Journal of Advances in Information Fusion, vol. 11, no. 2, pp.
%    138-156, Dec. 2016.
%[3] G. Nævdal, "On a generalization of the trigonometric moment
%    problem," Linear Algebra and Its Applications, vol. 258, pp. 1-18,
%    Jun. 1997.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(numComp))
        if(length(rhoList)==1)
            numComp=2;
        else
            numComp=5;
        end
    end
    
    if(any(abs(rhoList)>1))
        error('The magnitudes of the entries in rhoList must be <=1.')
    end
    
    switch(numComp)
        case 2%Section IIIA of [1] and 3.1.1 of [2].
            rho1=rhoList(1);
            mu=wrapRange(atan2(imag(rho1),real(rho1)),-pi,pi);
            phi=acos(abs(rho1));

            beta1=wrapRange(mu-phi,-pi,pi);
            beta2=wrapRange(mu+phi,-pi,pi);

            d=[beta1,beta2];
            w=[1/2;1/2];
        case 3%Section IIIA of [1] and Section 3.1.2 of [2].
            rho1=rhoList(1);
            mu=wrapRange(atan2(imag(rho1),real(rho1)),-pi,pi);
            phi=acos(3/2*abs(rho1)-1/2);
            beta1=wrapRange(mu-phi,-pi,pi);
            beta2=mu;
            beta3=wrapRange(mu+phi,-pi,pi);
            d=[beta1;beta2;beta3];
            w=[1/3;1/3;1/3];
        case 5%Section IIIB of [1] and Algorithm 1 in [2].
            if(nargin<3||isempty(lambda))
                lambda=0.8;
            end

            rho1=rhoList(1);
            rho2=rhoList(2);

            %We shall compute the determinant of 
            % X=[1,       rho1,     rho2;
            %   rho1',   1,        rho1;
            %   rho2',   rho1',    1];
            %but instead of using the actual rho2 value, since this
            %function does not match the phase of rho2, we use the rho2
            %value that has the same magnitude but the phase of rho2
            %optimized so that we get the maximum determinant possible.
            rho1R=real(rho1);
            rho1I=imag(rho1);
            detX=1-2*abs(rho1)^2-abs(rho2)^2+2*abs(rho2)*(rho1I^2+rho1R^2);

            if(detX<0)
               error('rhoList(2) has a value that does not lend itself to a solution.')
            end

            mu=wrapRange(atan2(imag(rho1),real(rho1)),-pi,pi);

            m1=abs(rho1);
            m2=abs(rho2);

            %Equation 8 in [1].
            w5Min=(4*m1^2-4*m1-m2+1)/(4*m1-m2-3);
            %Equation 9
            w5Max=(2*m1^2-m2-1)/(4*m1-m2-3);

            %Equation 10
            w5=w5Min+lambda*(w5Max-w5Min);

            %Defined in the text in Section IIIB in [1].
            w1=(1-w5)/4;

            %Equation 5
            c1=2/(1-w5)*(m1-w5);
            %Equation 6
            c2=1/(1-w5)*(m2-w5)+1;

            %Equation 7
            x2=(2*c1+sqrt(4*c1^2-8*(c1^2-c2)))/4;
            x1=c1-x2;

            %Deal with finite precision errors.
            x1=real(x1);
            x2=real(x2);

            %This error should have been avoided due to the check on the
            %derivative above.
            if(abs(x1)>1||abs(x2)>1)
                error('Unable to fit given moments. Please verify that rhoList(2) is representative of a real distribution.')
            end

            phi1=acos(x1);
            phi2=acos(x2);

            d=wrapRange(mu+[-phi1,-phi2,phi1,phi2,0],-pi,pi);
            w=[w1;w1;w1;w1;w5];
        otherwise
            error('Invalid number of components specified.')
    end
end

function vals=rand(N,d,w)
%%RAND Generate random samples from the wrapped Dirac distribution.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of  random variables.
%        d The 1XnumComp or numCompX1 vector of  real wrapped Dirac points.
%        w The numCompX1 weights associated with the wrapped Dirac points
%          sum(w)=1; w>=0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated wrapped Dirac random variables. These values range
%              from -pi to pi.
%
%Sampling a wrapped Dirac distribution is the same as sampling an empirical
%distirbution and wrapping the result to the range -pi to pi. Thus, this
%function just calls EmpiricalD.rand and wraps the result.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    vals=wrapRange(EmpiricalD.rand(N,d,w),-pi,pi);
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
