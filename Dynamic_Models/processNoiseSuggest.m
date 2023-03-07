function procNoiseParam=processNoiseSuggest(algorithm,maxVal,T,sigmaw2,manDur)
%%PROCESSNOISESUGGEST Use one of a number methods for choosing the scaling
%             parameter for the process noise covariance in a number of
%             continuous-time and discrete-time dynamic models. This only
%             provides the scaling parameter for one dimension of motion.
%             In a multidimensional system, if the measurement noise
%             variance varies between dimensions (such as tracking in range
%             and angle with a range-angle state), then one might use
%             different values in different dimensions. The method only
%             considers scalar measurements in the respective dimensions
%             and is thus only good for a rough approximation.
%
%INPUTS: algorithm A string specifying the algorithm and model used to
%           estimate a good process noise parameter. Possible values are:
%                  'PolyKal-ROT' Use a rule-of-thumb method for determining
%                                the q0 term for the QPolyKal function, for
%                                discrete-time models, and the DPoly
%                                function, for continuous-time models. This
%                                covers the continuous white noise
%                                acceleration (CWNA) model (order=1), the 
%                                discretized continuous white noise
%                                acceleration (DCWNA) model (order=1), the
%                                continuous Wiener process acceleration
%                                (CWPA) model  (order=2), and the 
%                                discretized continuous Wiener process
%                                acceleration (DCWPA) model (order=2),
%                                among others.
%           'PolyKalDirectDisc-ROT' Use a rule-of-thumb method for
%                                determining the sigmaV2 parameter for the
%                                QPolyKalDirectDisc function. This covers
%                                the discrete white noise acceleration
%                                (DWNA) model (order=1), among others.
%           'PolyKalDirectAlt-ROT' Use a rule-of-thumb method for
%                                determining the sigmaV2 parameter for the
%                                QPolyKalDirectAlt function. This covers
%                                the discrete Wiener process acceleration
%                                (DWPA) model (order=2), among others.
%                 'CWNA-OptMMSE' Use the method of Blair for determining
%                                the asymptotically optimal value of q0 in
%                                terms of MSE for the CWNA and DCWNA
%                                dynamic models for a maneuver having a
%                                fixed maximum acceleration.
%                 'DWNA-OptMMSE' Use the method of Blair for determining
%                                the asymptotically optimal value of
%                                sigmaV2 for the DWNA model for a maneuver
%                                having a fixed maximum acceleration.
%               'CWNA-ConstMeas' Use the method of Blair for determining
%                                the asymptotically optimal value of q0
%                                for the CWNA and DCWNA models such that
%                                the MSE of a Kalman filter under a maximum
%                                acceleration maneuver is not worse than
%                                the measurement accuracy.
%               'DWNA-ConstMeas' Use the method of Blair for determining
%                                the asymptotically optimal value of
%                                sigmaV2 for the DWNA model such that
%                                the MSE of a Kalman filter under a maximum
%                                acceleration maneuver is not worse than
%                                the measurement accuracy.
%    maxVal All of the methods require a maximum bound related to how the
%           target can maneuver. For algorithms CWNA-OptMMSE,
%           DWNA-OptMMSE, CWNA-ConstMeas, and DWNA-ConstMeas, maxVal is
%           the maximum acceleration of the target. For PolyKal-ROT and
%           PolyKalDirectDisc-ROT, and PolyKalDirectAlt-ROT maxVal is the
%           maximum value of a moment one order higher than the maximum
%           order of the dynamic model. For example, if order=1, then
%           maxVal is a maximum acceleration.  If order=2, then maxVal is
%           a maximum jerk.
%         T For all of the algorithms except PolyKalDirectDisc-ROT, this
%           parameter is required and is the typical time between
%           measurements of the target.
%   sigmaw2 For algorithms CWNA-OptMMSE, DWNA-OptMMSE, CWNA-ConstMeas, and
%           DWNA-ConstMeas, this parameter is required and is the variance
%           of the 1D measurement, which should be the lowest order moment
%           of the state. In other words, position.
%    manDur For algorithms CWNA-OptMMSE, DWNA-OptMMSE, CWNA-ConstMeas, and
%           DWNA-ConstMeas, this parameter is required and is the
%           maximum number of samples of duration T that a maximum
%           acceleration maneuver is expected to take. This can take
%           values, 3, 4, 5, 6 and infinity. The default if omitted or an
%           empty matrix is passed is Inf.
%
%OUTPUTS: procNoiseParam The value of q0 or sigmaV2 for the specified
%                        dynamic model, chosen according to the selected
%                        ad-hoc parameter.
%
%The output of this function goes into functions like QPolyKal. It is meant
%that one computes the process noise parameter once (for QPolyKal, this
%would correspond to the q0 parameter), and then one calls the other
%functions with different samples prediction intervals. This function
%should not be called to change the process noise parameter as the sampling
%period changes.
%
%All of the rule-of-thumb methods for choosing the process noise values are
%from Chapter 6.2 and 6.3 of [1] and modified slightly, as described in the
%comments for the implementation below.
%
%The solutions for the ConstMeas and MMSE methods are from [2], which is an
%extension of the work in [3] and [4].
%
%The rule-of thumb parameters might be suitable for other dynamic models.
%For example, when using the Singer model, given with aGaussMarkov and
%DPoly in continuous-time and with FGaussMarkov and QGaussMarkov in
%discrete-time, there does not appear to be a clear way to set the scaling
%value for the process noise in the literature. However, if tau=infinity in
%the model, then in continuous time, it reduces to the CWPA model. Thus,
%methods for choosing the process noise parameter in the CWPA and the DCWPA
%models might be good starting points for setting the process noise
%parameter in the continuous and discrete Singer models.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] W. D. Blair, "Design of nearly constant velocity filters for radar
%    tracking of maneuvering targets," in Proceedings of the IEEE Radar
%    Conference, Atlanta, GA, 7-11 May 2012, pp. 1008-1013.
%[3] W. D. Blair, "Design of nearly constant velocity track filter for
%    brief maneuvers," in Proceedings of the 14th International Conference
%    on Information Fusion, Chicago, IL, 5-8 Jul. 2011.
%[4] W. D. Blair, "Design of nearly constant velocity track filters for
%    tracking maneuvering targets," in Proceedings of the 11th
%    International Conference on Information Fusion, Cologne, Germany, 30
%    Jun. - 3 Jul. 2008.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(manDur))
    manDur=Inf;
end

switch(algorithm)
    case 'PolyKal-ROT'
        %The rule of thumb is generalized from the rules of thumb in
        %Chapters 6.2.2 and 6.2.3 of Bar-Shalom's book, when considering
        %the CWNA, DCWNA, CWPA and DCWPA. For the CWNA and DCWNA models, it
        %is suggested that sqrt(q0*T) be on the order of the changes in
        %velocity over a sampling interval. In other words, it is the
        %average acceleration. It might be more useful to express things in
        %terms of a maximum acceleration, though. In Chapter 6.3.2 for the
        %DWNA model, it is suggested that  0.5*aMax<=sigmaV<=aMax. Using a
        %similar logic, one might try  0.5*aMax<=sqrt(q0*T)<=aMax. Using
        %the same logic for the CWPA and DCWPA models of Chapter 6.2.3, one
        %might choose q0 such that  0.5*jerkMax<=sqrt(q0*T)<=jerkMax, where
        %jerkMax is the maximum possible jerk (derivative of acceleration).
        %This function generalizes the rule  to any order and chooses the
        %point midway in that specified range.
        
        procNoiseParam=(0.75*maxVal)^2/T;
        return;
    case 'PolyKalDirectDisc-ROT'
        %The rule-of-thumb is generalized from that given in Chapter 6.3.2
        %of Bar-Shalom's book. In the book, it is suggested that
        %0.5*aMax<=sigmaV<=aMax for the DWNA model (order=1). This function
        %generalizes it to any order and chooses the point midway in that
        %specified range.
        
        procNoiseParam=(0.75*maxVal)^2;
        return;
    case 'PolyKalDirectAlt-ROT'
        %The rule-of-thumb is generalized from that given in Chapter 6.3.3
        %of Bar-Shalom's book. In the book, it is suggested that
        %0.5*deltaAMax<=sigmaV<=deltaAMax for the DWPA model (order=2),
        %where deltaAMax is the maximum change in acceleration over an
        %interval of duration T. This can be related to a maximum jerk as
        %0.5*jerkMax*T<=sigmaV<=jerkMax*T. This function generalizes it to
        %any order and chooses the point midway in that range.
        
        procNoiseParam=(0.5*maxVal*T)^2;
        return;
    case 'CWNA-OptMMSE'
        optType=0;
        isDiscrete=false;
    case 'DWNA-OptMMSE'
        optType=0;
        isDiscrete=true;
    case 'CWNA-ConstMeas'
        optType=1;
        isDiscrete=false;
    case 'DWNA-ConstMeas'
        optType=1;
        isDiscrete=true;
    otherwise
        error('Unknown algorithm chosen')
end

%Equation 3, the deterministic maneuvering index.
GammaD=maxVal*T^2/sqrt(sigmaw2);

if(0.01>GammaD||GammaD>10)
   warning('The deterministic maneuvering index is outside of the tabulated range of values. The process noise parameter might be inaccurate.')
end
%Blair's papers do not say the base, but if one uses base 10, one gets the
%results in the plots. If one uses base e, the results do not match the
%plots in the paper.
logVal=log10(GammaD);

%optType=0 means MMSE (so use max)
%optType=1 means limit the maximum error to the measurement noise (so use
%          min)
switch(manDur)
    case Inf
        if(optType)
            curRow=2;
        else
            curRow=1;
        end
    case 3
        if(optType)
            curRow=4;
        else
            curRow=3;
        end
    case 4
        if(optType)
            curRow=6;
        else
            curRow=5;
        end
    case 5
    %If one wants a five-sample maneuver, then we will average the values
    %from the length four and the length-six maneuvers.
        if(optType)
            curRow=[6;8];
        else
            curRow=[5;7];
        end
    case 6
        if(optType)
            curRow=8;
        else
            curRow=7;
        end
    otherwise
        error('An untabulated value for the maximum maneuver duration has been given')
end

%Table I in Blair's paper
k1CoeffTable=[1.677, -0.726,  0.230, -0.012, 0.005,  0.000;%k^{max}_{\inf}
              0.872, -0.102, -0.019,  0.010, 0.006,  0.001;%k^{min}_{\inf}
              1.539, -0.189, -0.651,  0.187, 0.284,  0.059;%k^{max}_{3}
              0.707,  0.346, -0.318, -0.097, 0.088,  0.029;%k^{min}_{3}
              1.630, -0.473, -0.403,  0.272, 0.173,  0.023;%k^{max}_{4}
              0.803,  0.197, -0.385, -0.009, 0.138,  0.035;%k^{min}_{4}
              1.671, -0.697,  0.025,  0.250,-0.0324,-0.028;%k^{max}_{6}
              0.869,  0.008, -0.320,  0.088, 0.123,  0.023];%k^{min}_{6}

a0=k1CoeffTable(curRow,1);
a1=k1CoeffTable(curRow,2);
a2=k1CoeffTable(curRow,3);
a3=k1CoeffTable(curRow,4);
a4=k1CoeffTable(curRow,5);
a5=k1CoeffTable(curRow,6);

%Equations 6 and 7. The mean command is for averaging if one wanted a
%length 5 maneuver. Note that equation 6 has a repeated cubed term. That
%appears to be a mistake.
k1=mean(a0+logVal*(a1+logVal*(a2+logVal*(a3+logVal*(a4+a5*logVal)))));

if(isDiscrete)%Discrete-time
    procNoiseParam=(k1*maxVal)^2;
else%Continuous-time
    procNoiseParam=T*(k1*maxVal)^2;
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
