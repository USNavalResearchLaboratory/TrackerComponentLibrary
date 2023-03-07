function [xi,w]=fifthOrderCubPoints(numDim,algorithm,randomize)
%FIFTHORDERCUBPOINTS Generate fifth order cubature points for integration
%           involving a multidimensional Gaussian probability density
%           function (PDF).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated; numDim>=1.
%     algorithm An optional parameter specifying the algorithm to be used
%               to generate the points. Possible values are:
%               0 (The default if omitted or an empty matrix is passed
%                 and numDim~=2 and numDim~=1) Use algorithm E_n^{r^2} 5-3
%                 on page 317 of [1], requires 2^numDim+2*numDim points
%                 and numDim>=3.
%               1 Use formula I of [2], where some typos were corrected by
%                 referring back to [3]. Note that for numDim>7, some of
%                 the weights produced by the formula will be negative. It
%                 is required that numDim>=4.This formula uses fewer
%                 points than algorithm 1.
%               2 (The default if numDim==2) Use the cubature points first
%                 mentioned in [4], implemented as described in Section
%                 6.2.3 of [5]. The points are labeled as being fourth-
%                 order. However, due to their symmetry, they are actually
%                 fifth-order. The points are also given in [6].
%                 numDim>=2.
%               3 Formula E_n^{r^2} 5-2 on page 317 of [1],
%                 2*numDim^2+1 points. This is the same as the fifth-order
%                 formula derived in [7].
%               4 Formula E_n^{r^2} 5-4 on page 318 of [1], 2^(numDim+1)-1
%                 points.
%               5 Formula E_n^{r^2} 5-5 on page 318 of [1], choosing the
%                 upper signs, numDim*2^(numDim)+1 points.
%               6 Formula E_n^{r^2} 5-6 on page 318 of [1],
%                 2^numDim*(numDim+1) points, numDim>=5.
%               7 Formula E_2^{r^2} 5-1 on page 324 of [1], 7 points,
%                 numDim=2.
%               8 Formula E_2^{r^2} 5-2 on page 324 of [1], 9 points,
%                 numDim=2.
%               9 Formula E_3^{r^2} 5-1 on page 326 of [1], 13 points,
%                 numDim=3.
%              10 Formula E_3^{r^2} 5-2 on page 327 of [1], 14 points,
%                 numDim=3 (the first of two variants), corrected to take
%                 the fully symmetric permutations of the second vector.
%              11 Formula E_3^{r^2} 5-2 on page 327 of [1], 15 points,
%                 numDim=3 (the second of two variants), corrected to take
%                 the fully symmetric permutations of the second vector.
%              12 Formula E_3^{r^2} 5-3 on page 326 of [1], 21 points,
%                 numDim=3.
%              13 Formula E_4^{r^2} 5-1 on page 328 of [1], 25 points,
%                 numDim=4.
%              14 (The default if numDim==1) Call quadraturePoints1D(3), 3
%                 points, numDim=1.
%    randomize If this parameter is true, then the points will be
%              multiplied by a random orthonormal rotation matrix. This
%              does not change the moments up to the order of the points.
%              This randomization is done in [8] and [9] to lessen various
%              effects that arise when using points in the same orientation
%              repeatedly in tracking. The default if this parameter is
%              omitted or an empty matrix is passed is false. This
%              parameter is ignored if numDim==1.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This function returns fifth order cubature points and weights of the
%specified dimensions for Gaussian cubature integration over a normal
%distribution. The points and weights can be used to get exact solutions to
%multivariate integrals of the N(0,I) distribution times an polynomials up
%to the fifth order. When used with arbitrary functions, the points can be
%used to approximate the value of the integral.
%
%If f() if a function that returns the value for a function at a given
%point, then the integral of the N(0,I) times that function over the entire
%real space can be approximated using the points as
%
%intVal=0;
%for(i=1:numCubaturePoints)
%   intVal=intVal+w(i)*f(x(:,i));
%end
%
%If an integral over a function times N(mu,Sigma) is desired instead of
%over a function times N(0,I), the points should be modified as
%xi=transformCubPoints(xi,mu,chol(Sigma,'lower'))
%The weights remain unchanged.
%
%REFERENCES:
%[1] A. H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] J. Lu and D. L. Darmofal, "Higher-dimensional integration with
%    Gaussian weight for applications in probabilistic design," SIAM
%    Journal on Scientific Computing, vol. 26, no. 2, pp. 613-624, 2004.
%[3] I. P. Mysovskikh, "The approximation of multiple integrals by using
%    interpolatory cubature formulae," in Quantitative Approximation,
%    R. A. DeVore and K. Scherer, Eds. New York: Academic Press, 1980,
%    pp. 217-244.
%[4] S. J. Julier and J. K. Uhlmann, "A consistent, debiased method for
%    converting between polar and Cartesian coordinate systems," in
%    Proceedings of SPIE: Acquisition, Tracking, and Pointing XI, vol.
%    3086, Orlando, FL, 23-24 Oct. 1997, pp. 110-121.
%[5] U. N. Lerner, "Hybrid Bayesian networks for reasoning about complex
%    systems," Ph.D. dissertation, Stanford University, Stanford, CA,
%    Oct. 2002.
%[6] Y. Wu, D. Hu, M. Wu, and X. Hu, "A numerical-integration perspective
%    on Gaussian filters," IEEE Transactions on Signal Processing, vol. 54,
%    no. 8, pp. 2910-2921, Aug. 2006.
%[7] B. Jia, M. Xin, and Y. Cheng, "High-degree cubature Kalman filter,"
%    Automatica, vol. 49, no. 2, pp. 510-518, Feb. 2013.
%[8] O. Straka, D. Duník, and M. Simandl, "Randomized unscented Kalman
%    filter in tracking," in Proceedings of the 15th International
%    Conference on Information Fusion, Singapore, 9-12 Jul. 2012, pp.
%    503-510.
%[9] J. Duník, O. Straka, and M. Simandl, "The development of a randomised
%    unscented Kalman filter," in Proceedings of the 18th World Congress,
%    The International Federation of Automatic Control, Milan, Italy, 28
%    Aug. - 2 Sep. 2011, pp. 8-13.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(algorithm))
        if(numDim==2)
            algorithm=2;
        elseif(numDim==1)
            algorithm=14;
        else
            algorithm=0;
        end
    end
    
    if(nargin<3||isempty(randomize))
        randomize=false;
    end

    n=numDim;
    switch(algorithm)
        case 0%Algorithm 5-3 on page 317 of [1].
            if(n<3)
                error('This algorithm requires numDim>=3.')
            end
            
            numPoints=2^n+2*n;

            w=zeros(numPoints,1);
            w(1:(2*n))=4/(n+2)^2;
            w((2*n+1):numPoints)=(n-2)^2/(2^n*(n+2)^2);

            %Allocate space
            xi=zeros(n,2*n);

            %First, go through all of the positions of r:
            r=sqrt((n+2)/2);

            %First, the positive values
            for k=1:n
                xi(k,k)=r; 
            end

            %Now, the negative values
            for k=(n+1):(2*n)
                xi(k-n,k)=-r;
            end

            %Now, go through all of the positions of s.
            s=sqrt((n+2)/(n-2));
            xi=[xi, PMCombos(ones(n,1)*s,0)];
        case 1
            if(n<4)
                error('This algorithm requires numDim>=4.')
            end
            numPoints=n^2+3*n+3; 
            w=zeros(numPoints,1);
            xi=zeros(n,numPoints);

            xCur=1;
            w(xCur)=2/(n+2);%The first point is the origin.

            %The next 2*(n+1) points come from calculating a-values, which
            %must be saved to calculate b-values.
            ub=2*(n+1)+xCur;
            xCur=xCur+1;
            w(xCur:ub)=n^2*(7-n)/(2*(n+1)^2*(n+2)^2);

            ar=zeros(n,n+1);
            temp=sqrt(n+2);
            for r=1:(n+1)
                for k=1:(r-1)
                    ar(k,r)=-sqrt((n+1)/(n*(n-k+2)*(n-k+1)));
                end
                %for k=r;
                %This had a TYPO in Lu and Darmofal's paper that was 
                %correcting using Mysovskikh paper (which itself is missing
                %a closed parenthesis in one part).
                if(r<n+1)
                    ar(r,r)=sqrt((n+1)*(n-r+1)/(n*(n-r+2)));
                end
                %for all higher k, it is zero.

                xi(:,xCur)=temp*ar(:,r);
                xi(:,xCur+1)=-temp*ar(:,r);
                xCur=xCur+2;
            end

            %Now, the final n*(n+1) points come from combining all pairs of
            %points.
            w(xCur:end)=2*(n-1)^2/((n+1)^2*(n+2)^2);

            temp=sqrt(0.5*(n*(n+2)/(n-1)));
            for l=1:n
                for k=(l+1):(n+1)
                    b=temp*(ar(:,l)+ar(:,k));
                    xi(:,xCur)=b;
                    xi(:,xCur+1)=-b;
                    xCur=xCur+2;
                end
            end
        case 2
            if(n<2)
                error('This algorithm requires numDim>=2.')
            end
        %The numbers of different types of points.
            numType1=1;
            numType2=2*n;
            numType3=2*n*(n-1);
            numPoints=numType1+numType2+numType3;

        %First, calculate the weights.
            w=zeros(numPoints,1);
            w(1)=(1/18)*(18-7*n+n^2);
            w((numType1+1):(numType1+numType2))=(4-n)/18;
            w((numType1+numType2+1):numPoints)=1/36;

        %Now, calculate the sigma points.
            xi=zeros(n,numPoints);

            %The first point is a type-1 point and is left at zero.

            %The type two points come in +/- pairs on the coordinate axes.The
            %value on the coordinate axes is
            u=sqrt(3);
            count=2;
            for idx=1:n
                xi(idx,count)=u;
                xi(idx,count+1)=-u;

                count=count+2;
            end
    
            %The type 3 points are formed by going through all possible
            %combinations of the data such that (if n>2), one element is
            %zero and all of the other elements are either u or -u. The
            %points cover all such combinations.

            %We need to go through all possible combinations of choosing
            %two spots in a vector of length n to place the values. the
            %values can be +/-u.

            combos=nchoosek(1:n,2)';
            numCombos=size(combos,2);

            for curCombo=1:numCombos
                curSel=combos(:,curCombo);
                xi(curSel,count)=[u;u];
                xi(curSel,count+1)=[u;-u];
                xi(curSel,count+2)=[-u;u];
                xi(curSel,count+3)=[-u;-u];
                count=count+4;
            end
        case 3%Formula E_n^{r^2} 5-2 on page 317 of [1], 2*numDim^2+1
              %points.

            A=2/(n+2);
            B=(4-n)/(2*(n+2)^2);
            C=1/(n+2)^2;
            
            r=sqrt((n+2)/2);
            s=sqrt((n+2)/4);
            
            %The sqrt(2) coefficient makes it for a normal 0-I
            %distribution.
            xi=sqrt(2)*[zeros(numDim,1),fullSymPerms([r;zeros(numDim-1,1)]),fullSymPerms([s;s;zeros(numDim-2,1)])];
            w=[A;B*ones(2*n,1);C*ones(2*n*(n-1),1)];
        case 4%Formula E_n^{r^2} 5-4 on page 318 of [1], 2^(numDim+1)-1
              %points.
            numPoints=2^(n+1)-1;
            xi=zeros(numDim,numPoints);
            w=zeros(numPoints,1);
            
            s=1/sqrt(2);
            
            %Zeroth-order term has xi all zeros and this for w
            w(1)=2/(n+2);
            curStartPoint=2;
            
            %Build the first order term.
            theVec=s*ones(numDim,1);
            for k=1:n
                theVec(k)=sqrt((k+2)/2);
                xiCur=PMCombos(theVec(k:end));
                num2Add=size(xiCur,2);
                
                xi(k:end,curStartPoint:(curStartPoint+num2Add-1))=xiCur;
                w(curStartPoint:(curStartPoint+num2Add-1))=2^(k-n)/((k+1)*(k+2));

                curStartPoint=curStartPoint+num2Add;
            end
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*xi;
        case 5%Formula E_n^{r^2} 5-5 on page 318 of [1], choosing the upper
              %signs, numDim*2^(numDim)+1 points
            r=sqrt((n+2+(n-1)*sqrt(2*(n+2)))/(2*n));
            s=sqrt((n+2-sqrt(2*(n+2)))/(2*n));
            
            A=(2/(n+2));
            B=1/(2^n*(n+2));
            
            numPoints=n*2^n+1;
            xi=zeros(n,numPoints);
            w=zeros(numPoints,1);
            
            %The first point is all zeros.
            w(1)=A;
            
            %Filling in othe other points and scaling by sqrt(2) to make it
            %for the normal 0-I distribution.
            xi(:,2:end)=sqrt(2)*fullSymPerms([r;s*ones(n-1,1)]);
            w(2:end)=B;
        case 6%Formula E_n^{r^2} 5-6 on page 318 of [1],
              %2^numDim*(numDim+1) points, numDim>=5.
            if(n<5)
                error('This algorithm requires numDim>=5.')
            end
            
            r=sqrt((n-sqrt(2)+(n-1)*sqrt(2*(n+1)))/(2*n));
            s=sqrt((n-sqrt(2)-sqrt(2*(n+1)))/(2*n));
            t=sqrt((1+sqrt(2))/2);

            numPoints=2^n*(n+1);
            
            A=1/numPoints;
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*[fullSymPerms([r;s*ones(n-1,1)]),PMCombos(t*ones(n,1))];
            w=A*ones(numPoints,1);
        case 7%Formula E_2^{r^2} 5-1 on pg. 324 of [1] 7 points, numDim=2.
            if(n~=2)
                error('This algorithm requires numDim=2.')
            end
            
            r=sqrt(2);
            s=sqrt(2)/2;
            t=sqrt(6)/2;
            A=1/2;
            B=1/12;
            
            w=[A;B*ones(6,1)];
            xi=zeros(2,7);
            %xi(:,1)is zero.
            
            xi(:,2)=[r;0];
            xi(:,3)=[-r;0];
            xi(:,4)=[s;t];
            xi(:,5)=[s;-t];
            xi(:,6)=[-s;t];
            xi(:,7)=[-s;-t];
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*xi;
        case 8%Formula E_2^{r^2} 5-2 on pg. 324 of [1], 9 points, numDim=2.
            if(n~=2)
                error('This algorithm requires numDim=2.')
            end
            
            A=4/9;
            B=1/9;
            C=1/36;
            
            r=sqrt(3/2);
            
            w=[A;B*ones(4,1);C*ones(4,1)];
            xi=zeros(2,9);
            %xi(:,1)is zero.
            
            xi(:,2:5)=fullSymPerms([r;0]);
            xi(:,6:9)=PMCombos([r;r]);

            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*xi;
        case 9%Formula E_3^{r^2} 5-1 on page 326 of [1], 13 points,
              %numDim=3;
            if(n~=3)
                error('This algorithm requires numDim=3.')
            end
            
            A=2/5;
            B=1/20;
            r=sqrt((5+sqrt(5))/4);
            s=sqrt((5-sqrt(5))/4);
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(3,1),PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r])];
            w=[A;B*ones(12,1)];
        case 10%Formula E_3^{r^2} 5-2 on page 327 of [1], 14 points,
               %numDim=3 (the first of two variants), corrected to take the
               %fully symmetric permutations of the second vector.
            if(n~=3)
                error('This algorithm requires numDim=3.')
            end
            
            r=sqrt(5/4);
            s=sqrt(5/2);
            B=4/25;
            C=1/200;
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*[fullSymPerms([r;0;0]),PMCombos([s;s;s])];
            w=[B*ones(6,1);C*ones(8,1)];
        case 11%Formula E_3^{r^2} 5-2 on page 327 of [1], 15 points,
               %numDim=3 (the second of two variants), corrected to take the
               %fully symmetric permutations of the second vector.
            if(n~=3)
                error('This algorithm requires numDim=3.')
            end
            
            r=sqrt(5/2);
            s=sqrt(5/6);
            A=2/5;
            B=1/25;
            C=9/200;
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*[[0;0;0],fullSymPerms([r;0;0]),PMCombos([s;s;s])];
            w=[A;B*ones(6,1);C*ones(8,1)];
        case 12%Formula E_3^{r^2} 5-3 on page 327 of [1], 21 points,
               %numDim=3.
            if(n~=3)
                error('This algorithm requires numDim=3.')
            end
            
            r=sqrt((15-5*sqrt(5))/12);
            s=sqrt((15+5*sqrt(5))/12);
            t=sqrt(5/6);
            A=2/5;
            B=3/100;
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(3,1),PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([t;t;t])];
            w=[A;B*ones(20,1)];
        case 13%Formula E_4^{r^2} 5-1 on page 328 of [1], 25 points,
               %numDim=4.
            if(n~=4)
                error('This algorithm requires numDim=4.')
            end
            
            A=1/3;
            B=1/36;
            r=sqrt(3/2);
            
            %Scaling by sqrt(2) to make it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(4,1),fullSymPerms([r;r;0;0])];
            w=[A;B*ones(24,1)];
        case 14%Use quadraturePoints1D(3) for 1D points
            if(n~=1)
                error('This algorithm requires numDim=1.')
            end
            [xi,w]=quadraturePoints1D(3);
        otherwise
            error('Unknown cubature algorithm selected')
    end
    
    if(numDim>1&&randomize)
        R=randOrthoMat(numDim);
        
        numPoints=length(w);
        for curPoint=1:numPoints
            xi(:,curPoint)=R*xi(:,curPoint);
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
