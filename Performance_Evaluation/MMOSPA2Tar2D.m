function MMOSPA=MMOSPA2Tar2D(particles)
%%MMOSPA2TAR2D Given a set of uniformly weighted particles representing the
%       positions of two targets in 2D space, this function find the exact
%       minimum mean optimal subpattern assignment (MMOSPA) estimate of the
%       location of each of the targets using the efficient polynomial
%       algorithm of [1].
%
%INPUTS: particles A 4XnumParticles set of particles representing two 2D
%           sets of target positions where each particle has the first
%           targets components followed by those of the second target. This
%           is:
%           [x1;y1;x2;y2]
%
%OUTPUTS: MMOSPA The 4X1 MMOSPA estimate where MMOSPA(1:2) are the
%                components of the first target and MMOSPA(3:4) are the
%                components of the second target.
%
%This implements algorithm 4 in [1] assuming that all of the weights are
%uniform.
%
%EXAMPLE:
% numParticles=1000;
% %The distribution of the targets is modeled here as a Gaussian mixture with
% %two components in 2D. The mean is 4X2 because the dimensionality of the
% %stacked state is 4D and there are two components in the mixture used here.
% mu=[2,4;
%     2,0;
%     4,2;
%     0,2];
% %A 4X4X2 matrix of joint covariance matrices for the two mixture
% %components.
% P=cat(3,0.5*eye(4),0.3*eye(4));
% w=[1/3,2/3];%Mixture weights.
% particles=GaussianMixtureD.rand(numParticles,w,mu,P);
% 
% %Compute the estimates.
% MMSE_estimate=mean(particles,2);
% MMOSPA_estimate=MMOSPA2Tar2D(particles);
% %Plotting
% figure()
% axis equal
% hold on
% %Plot particles
% plot(particles(1,:), particles(2,:),'o','MarkerEdgeColor','black', 'MarkerFaceColor','magenta','MarkerSize',4);
% plot(particles(3,:), particles(4,:),'o','MarkerEdgeColor','black', 'MarkerFaceColor','blue', 'MarkerSize',4 );
% %Plot Estimates
% plot(MMSE_estimate([1;3]),MMSE_estimate([2;4]),'x','MarkerSize',12,'LineWidth',3,'Color','red');
% plot(MMOSPA_estimate([1;3]),MMOSPA_estimate([2;4]),'o','MarkerSize',12,'LineWidth',2,'Color','black','MarkerFaceColor','green');
% legend('Particles for Target 1', 'Particles for Target 2', 'MMSE Estimate','MMOSPA Estimate');
% xlabel('x');
% ylabel('y');
%
%REFERENCES:
%[1] M. Baum, P. Willett, and U. D. Hanebeck, "Polynomial-time algorithms
%    for the exact MMOSPA estimate of a multi-object probability density
%    represented by particles," IEEE Transactions on Signal Processing,
%    vol. 63, no. 10, 15 May 2015, pp. 2476-2484.
%
%September 2023 Marcus Baum, University of Gottingen, Gottingen, Germany.
%Modified (September 2023) for the Tracker Component Library by David F.
%Crouse, Naval Research Laboratory, Washington D.C., 
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(particles,2);
numDim=size(particles,1);

if(numDim~=4)
    error('This function only works with two 2D states.')
end

% Sort particles according to their orientation 
d=(particles(1:2,:)-particles(3:4,:));
dir=atan(d(2,:)./d(1,:));
[~,perm_dir]=sort(dir,2);
particles=particles(:,perm_dir);
% Initial permutation of particles
for j=1:n
   if([-1,0]*(particles(1:2,j)-particles(3:4,j))< 0)
      particles(:,j)=particles([3,4,1,2],j);
   end 
end
%Find Minimum MOSPA
MOSPA_temp=mean(particles,2);
reward_temp=MOSPA_temp'*MOSPA_temp;
MMOSPA=MOSPA_temp;
max_reward=reward_temp; 
for j=1:n
    switched_particle=particles([3,4,1,2],j)/n;
    temp=MOSPA_temp-particles(:,j)/n; 
    MOSPA_temp=temp+switched_particle;
    reward_temp=temp'*temp+2*temp'*switched_particle+switched_particle'*switched_particle;
    if(reward_temp>max_reward)
       max_reward=reward_temp;
       MMOSPA=MOSPA_temp;   
    end
end
end

%The MIT License:
%Copyright 2023, Marcus Baum (University of Gottingen)
% 
%Permission is hereby granted, free of charge, to any person obtaining a
%copy of this software and associated documentation files (the "Software"),
%to deal in the Software without restriction, including without limitation
%the rights to use, copy, modify, merge, publish, distribute, sublicense,
%and/or sell copies of the Software, and to permit persons to whom the
%Software is furnished to do so, subject to the following conditions:
% 
%The above copyright notice and this permission notice shall be included in
%all copies or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
%THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
%DEALINGS IN THE SOFTWARE.
