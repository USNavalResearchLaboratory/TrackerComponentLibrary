function beta=calc2DAssignmentProbsApprox(A,approxType,diagAugment,delta)
%%CALC2DASSIGNMENTPROBSAPPROX Given a matrix of all-positive likelihoods or
%                likelihood ratios, determine the probability that each row
%                (target) is assigned to each column (measurement or missed
%                detection hypothesis). Multiple approximation algorithms
%                are available. This function is useful when the function
%                calc2DAssignmentProbs is too slow. Missed detection
%                likelihoods are required (but can be zero) and are
%                on to the main diagonal of a diagonal matrix appended to
%                the end of the likelihood matrix.
%
%INPUTS:    A  A matrix of positive likelihoods or likelihood ratios (NOT
%              log-likelihood ratios). If diagAugment=true (required for
%              all approximations except approxTyp=0), then A is a
%              numTar X  (numMeas+numTar) matrix of all-positive
%              likelihoods or likelihood ratios for assigning the target
%              specified by the row to the measurement specified by the
%              column. Columns > numMeas hold missed-detection likelihoods/
%              likelihood ratios. Thus, off-diagonal terms for columns >
%              numMeas should be set to 0 and the diagonal terms set to the
%              costs of a missed detection for each given target.
%   approxType An optional parameter specifying the type of approximation
%              to use. The default if omitted is 0. Possible values are
%              0 The sum-product belief propagation algorithm of [6].
%              1 The Cheap JPDA of [1].
%              2 The modified cheap JPDAF of [2].
%              3 The algorithm of [3].
%              4 The matrix permanent approximation algorithm of [4].
%  diagAugment A boolean variable indicating whether the probabilties
%              should be found for a general assignment problem or assuming
%              A ends with a diagonal submatrix of missed detection
%              probabilities. The default if omitted is false (the general
%              problem). Setting diagAugment to true changes the shape of
%              the output. The default if omitted is true. See the
%              description of the output beta for how this affects the
%              output.
%        delta The optional convergence parameter used by the approxType=0.
%              If omitted or an empty matrix is passed, the default of 1e-4
%              is used. This is roughly proportional to the number of
%              digits of precision for convergence (not accuracy).
%
%OUTPUTS: beta  beta is a numTar X (numMeas+1) (if diagAugment=true) or a
%               numTar X (numMeas) (if diagArgument=false) matrix of
%               probabilities of assigning the target given by the row to
%               the measurement given by the column. The numMeas+1 column,
%               if present, is a set of missed detection probabilities.
%               Thus, in this case, beta has fewer columns than A, because
%               the final numRow columns are collapsed into one column on
%               return.
%
%Approximations to target-measurement association hypotheses are necessary
%when a large number of targets/ measurements are present as the complexity
%scales exponentially with the number of targets. The references below are
%the algorithms; reference [5] compares the algorithms in [1]-[4].
%
%REFERENCES:
%[1] R. J. Fitzgerald, "Development of practical PDA logic for multitarget
%    tracking by microprocessor," in Multitarget-Multisensor Tracking:
%    Advanced Applications, Y. Bar-Shalom, Ed. Norwood, MA.: Artech House,
%    1996, pp. 1-23.
%[2] P. Quan, Z. Hongcai, W. Peide and H. Zhou, "Development of new
%    practical multitarget data association algorithm," Proceedings of the
%    First IEEE Conference on Control Applications, vol. 2, Dayton, OH,
%    Sep. 1992, pp. 1065-1066.
%[3] B. Bakhtiar and H. Alavi, "Efficient algorithm for computing data
%    association probabilities for multitarget tracking," Proceedings of
%    SPIE: Automatic Object Recogniation IV, vol. 2756, Apr. 1996, pp.
%    130-140.
%[4] J. K. Uhlmann, "Matrix permanent inequalities for approximating joint
%    assignment matrices in tracking systems," Journal of the Franklin
%    Institute, vol. 341, no. 7, pp. 569-593, Nov. 2004.
%[5] K. Romeo, D. F. Crouse, Y. Bar-Shalom, and P. Willett, "The JPDAF in
%    practical systems: Approximations," in Proceedings of SPIE:
%    Signal and Data Processing of Small Targets Conference, vol. 7698,
%    Apr. 2010.
%[6] J. L. Williams and R. A. Lau, "Approximate evaluation of marginal
%    association probabilities with belief propagation," IEEE Transactions
%    on Aerospace and Electronic Systems, vol. 50, no. 4, pp. 2942-2959,
%    Oct. 2014.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(approxType))
    approxType=0;
end

if(nargin<3||isempty(diagAugment))
    diagAugment=true;
end

if(nargin<4||isempty(delta))
    delta=1e-3;
end

numTar=size(A,1);

%The algorithms other than the belief propagation method assume that the
%matrix is augmented.
if(diagAugment==false && approxType~=0)
    error('Only approximation 0 is suitable for use with problems involving no missed detection hypotheses')
end

numMeas=size(A,2)-numTar;
beta=zeros(numTar,numMeas+1);

switch(approxType)
    case 0%Sum-product belief propagation method
        beta=sumProdAssocProbApprox(A,diagAugment,delta);
    case 1%Cheap JPDAF
        for curTar=1:numTar
            Str=sum(A(curTar,1:numMeas));
            %Missed detection likelihood for this target
            B=A(curTar,numMeas+curTar);
            for curMeas=1:numMeas 
                Sji=sum(A(1:numTar,curMeas));

                Gtj=A(curTar,curMeas);

                beta(curTar,curMeas)=Gtj/(Sji+Str-Gtj+B);
            end

            %The missed detection hypothesis
            beta(curTar,numMeas+1)=1-sum(beta(curTar,1:numMeas));
        end
    case 2%Quan, Hongcai, Peide and Zhou's Algorithm
        %First, find the r's.
        r=zeros(numTar,numMeas);
        for curTar=1:numTar
            Str=sum(A(curTar,1:numMeas));
            for curMeas=1:numMeas
                r(curTar,curMeas)=A(curTar,curMeas)/Str;
            end
        end

        %Next, compute the approximate probabilities.
        for curTar=1:numTar
            sum1=sum(A(curTar,1:numMeas));
            for curMeas=1:numMeas
                sum2=sum(r(:,curMeas).*A(:,curMeas));
                %Missed detection likelihood for this target
                B=A(curTar,numMeas+curTar);
                beta(curTar,curMeas)=A(curTar,curMeas)/(sum1+sum2-r(curTar,curMeas)*A(curTar,curMeas)+B); 
            end
            
            %The missed detection hypothesis
            beta(curTar,numMeas+1)=1-sum(beta(curTar,1:numMeas));
        end
    case 3%The Bakhtiar and Alavi Algorithm
        beta=zeros(numTar,numMeas+1);
        
        %It is simpler just to put the missed detection hypotheses as an
        %extra column at the end.
        GMissedDet=diag(A(:,(numMeas+1):end));
        G=[A(:,1:numMeas),GMissedDet];
        for curTar=1:numTar
           tarSel=[1:(curTar-1),(curTar+1):numTar];

           for curMeas=1:numMeas
               %The +1 selects the missed detection likelihood
               measSel=[1:(curMeas-1),(curMeas+1):(numMeas+1)];
               
               beta(curTar,curMeas)=G(curTar,curMeas)*prod(sum(G(tarSel,measSel),2));
           end
           %The missed detection hypothesis
           beta(curTar,numMeas+1)=G(curTar,curMeas+1)*prod(sum(G(tarSel,:),2));
        end
        
        %Normalize the probabilities.
        beta=bsxfun(@rdivide,beta,sum(beta,2));
    case 4%Uhlmann's algorithm using approximate matrix permanents.
        numTar=size(A,1);
        numCol=size(A,2);
        numMeas=size(A,2)-numTar;

        boolRowsSkip=false(numTar,1);
        boolColsSkip=false(numCol,1);

        beta=zeros(numTar,numMeas+1);
        for curTar=1:numTar
            boolRowsSkip(curTar)=true;

             for curMeas=1:numMeas
                 ati=A(curTar,curMeas);
                 %The A matrix removing the row and column corresponding to the
                 %selected association.
                 boolColsSkip(curMeas)=true;
                 beta(curTar,curMeas)=ati*permBound(A(~boolRowsSkip,~boolColsSkip));

                 boolColsSkip(curMeas)=false;
             end
             %The missed detection hypothesis
             curMeas=numMeas+curTar;
             ati=A(curTar,curMeas);

             boolColsSkip(curMeas)=true;
             beta(curTar,numMeas+1)=ati*permBound(A(~boolRowsSkip,~boolColsSkip));

             boolColsSkip(curMeas)=false;
             boolRowsSkip(curTar)=false;
        end

        %Normalize the probabilities.
        beta=bsxfun(@rdivide,beta,sum(beta,2));
    otherwise
        error('Invalid Approximation Specified')
end
end

function beta=sumProdAssocProbApprox(A,diagAugment,delta)
%%SUMPRODASSOCPROBAPPROX Given a matrix of all-positive likelihoods or
%                likelihood ratios, determine the probability that each row
%                (target) is assigned to each column (measurement or missed
%                detection hypothesis) using an approximation based on the
%                sum-product algorithm. This function approximates the
%                function calc2DAssignmentProbs but is usually much faster.
%
%INPUTS:    A  A matrix of positive likelihoods or likelihood ratios (NOT
%              log-likelihood ratios). If diagAugment=true, then A is a
%              numTar X  (numMeas+numTar) matrix of all-positive
%              likelihoods or likelihood ratios for assigning the target
%              specified by the row to the measurement specified by the
%              column. Columns > numMeas hold missed-detection likelihoods/
%              likelihood ratios. Thus, off-diagonal terms for columns >
%              numMeas should be set to 0 and the diagonal terms set to the
%              costs of a missed detection for each given target.
%  diagAugment A boolean variable indicating whether the probabilties
%              should be found for a general assignment problem or assuming
%              A ends with a diagonal submatrix of missed detection
%              probabilities. The default if omitted is false (the general
%              problem). Setting diagAugment to true changes the shape of
%              the output. The default if omitted is false. See the
%              description of the output beta for how this affects the
%              output.
%        delta The optional convergence parameter used by the algorithm as
%              in [2]. If omitted or an empty matrix is passed, the default
%              of 1e-4 is used. This is rougly proportional to the number
%              of digits of precision for convergence (not accuracy).
%
%OUTPUTS:  beta If diagAugment is omitted or false, then beta has the same
%               dimensionality as A and hold the probability of assigning
%               each row to each column in the traditional 2D assignment
%               problem (Each row must be assigned to at least one column,
%               each column can be assigned to at most one row, or vice
%               versa if there are more rows than columns. If diagAugment
%               is true, then beta is a numTar X (numMeas+1) matrix of
%               probabilities of assigning the target given by the row to
%               the measurement given by the column. The final column is a
%               set of missed detection probabilities. Thus, in this case,
%               beta has fewer columns than A, because the final numRow
%               columns are collapsed into one column on return.
%
%This implements the algorithm described in [2], which is also presented in
%a slightly different manner in [1]. That algorithm is appropriate for the
%case where diagAugment=true. For the case where there are no missed
%diagAugment=false, the 1 in Equation 20 in [2] is removed.
%
%When diagAugment=true, the algorithm is run as in [2] and it is assumed
%that none of the missed detection likelihood is zero. This assumption
%is necessary, because the algorithm in [2] requires that the cost of the
%missed detection for each target (row) be 1, so the likelihoods in each
%row are divided by their relevant missed detection likelihood.
%
%REFERENCES:
%[1] J. L. Williams and R. A. Lau, "Convergence of loopy belief propagation
%    for data association," in Proceedings of the Sixth International
%    Conference on Intelligent Sensors, Sensor Networks and Information
%    Processing, Brisbane, Australia, 7-10 Dec. 2010, pp. 175-180.
%[2] J. L. Williams and R. A. Lau, "Approximate evaluation of marginal
%    association probabilities with belief propagation," IEEE Transactions
%    on Aerospace and Electronic Systems, vol. 50, no. 4, pp. 2942-2959,
%    Oct. 2014.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(diagAugment))
   diagAugment=false; 
end

if(nargin<3||isempty(delta))
   delta=1e-3;%Convergence criterion
end

n=size(A,1);
m=size(A,2);

if(isempty(A))
    beta=[];
    return;
end

if(diagAugment)
    %We will divide through the missed detection probabilities so that the
    %likelihood of a missed detection is 1.
    numMeas=m-n;
    
    %If the only hypothesis is the missed detection hypothesis because
    %there are no measurements.
    if(numMeas==0)
        beta=ones(n,1);
        return;
    end
    
    missedDetectProbs=diag(A(:,(numMeas+1):m));
    A=A(:,1:numMeas);
    A=bsxfun(@times,A,1./missedDetectProbs);
    m=numMeas;
end

Psi=A;

%Number of iterations to perform each step before checking for convergence.
N=100;%In [2], between 5 and 20 was recommended.

WStar=max(sum(Psi,2));%For determining convergence.

mu=zeros(n,m);%Target to measurement messages
v=100*ones(m,n);%Measurement to target messages, initialized to one.

if(diagAugment)
%If an implied missed detection hypothesis with cost 1 for all targets
%should be added.
    while(1)
        %Do a block of N iterations before checking convergence
        for k=1:N
            %Compute target to measurement messages
            for i=1:n
                prodTerms=Psi(i,:).*v(:,i)';

                s=1+sum(prodTerms);
                mu(i,:)=Psi(i,:)./(s-prodTerms);
            end

            %If it is the last iteration where a convergence check must be
            %performed
            if(k==N)
               vTilde=v; 
            end

            %Compute measurement to target messages
            for j=1:m
                s=1+sum(mu(:,j));
                v(j,:)=1./(s-mu(:,j));
            end
        end
        %Check for convergence
        d=max(max(abs(log(v./vTilde))));
        alpha=log((1+WStar*d)/(1+WStar))/log(d);

        if(alpha*d/(1-alpha)<=(1/2)*log(1+delta))
           break; 
        end
    end

    %Compute the target-measurement association probabilities.
    beta=zeros(n,m+1);
    for i=1:n
        prodTerms=Psi(i,:).*v(:,i)';
        s=1+sum(prodTerms);

        beta(i,1:m)=prodTerms/s;
        beta(i,m+1)=1/s;
    end
else
%If all hypotheses are represented and there is no missed detection
%hypothesis.
    while(1)
        %Do a block of N iterations before checking convergence
        for k=1:N
            %Compute target to measurement messages
            for i=1:n
                prodTerms=Psi(i,:).*v(:,i)';

                s=sum(prodTerms);
                mu(i,:)=Psi(i,:)./(s-prodTerms(i));
                
                %Deal with some finite precision issues
                sel=~isfinite(mu(i,:));
                mu(i,sel)=0;
            end

            %If it is the last iteration where a convergence check must be
            %performed
            if(k==N)
               vTilde=v; 
            end

            %Compute measurement to target messages
            for j=1:m
                s=1+sum(mu(:,j));
                v(j,:)=1./(s-mu(:,j));
            end
        end
        %Check for convergence
        d=max(max(abs(log(v./vTilde))));
        alpha=log((1+WStar*d)/(1+WStar))/log(d);

        if(alpha*d/(1-alpha)<=(1/2)*log(1+delta))
           break; 
        end
    end

    %Compute the target-measurement association probabilities.
    beta=zeros(n,m);
    for i=1:n
        prodTerms=Psi(i,:).*v(:,i)';
        s=sum(prodTerms);

        beta(i,1:m)=prodTerms/s;        
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
