function [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(Deltat,N,L,M,m,mu,gamma,W,fitntype,s,K,deltamax,C,tvec,qvec,tunnel)
% For a number of iterated scenarios, evaluate statistics of the waiting time T until m out of M 
% requlatory sequences of length L acquire a binding site of length W, so that the expression of the
% corresponding m genes is changed. The computed statistics, for all iterated scenarios, are the
% expected value the standard deviation, the density function, and quantile function of T.
%
% More extensive documentation is available in waitingtime_regseq_whaleevolution.m, where the waiting 
% time is computed for each specific iteration. The number of iterations (nriter) is computed from the 
% dimensionality of those input parameters that can be varied between iterations.

% Input parameters
% Deltat   - 1xnriter vector of generation time (typically in years)
% N        - 1xnriter vector of popoulation sizes
% L        - 1xnriter vector of lengths of the regulatory regions of all genes
% M        - 1xnriter vector of candidate genes to be part of the final target  
% m        - 1xnriter vector of number of genes that make up a target
% mu       - 1xnriter vector of mutation rates per generation and locus
% gamma    - 1xnriter vector of backward mutation probabilities
% W        - 1xnriter vector of lengths of binding sites
% fitntype - Specifies fitness function, either: 
%            'Mult' (multiplicative), 
%            'Final' (final target has different fitness)
%            'Stair' (intermediate type and final targets have different fitness)
%            'General' (general type of fitness funcrtion) 
% s        - Specifies fitness parameters. The input format is either:
%            fitntype = 'Mult' or 'Final': 1xnriter of selection parameters
%            fitntype = 'Stair': nriter x 2 matrix of selection parameters
%            fitntype = 'General': nriter x (max(M)+1) matrix of selection coefficients, where 
%                        s(i,j)=1 if j > M(i)+1.
% K        - 1xnriter vector of possible binding sites per gene
% deltamax - 1xnriter vector of maximum number of mismatches
% C        - 1xnriter vector of number of distance intervals per gene
% tvec     - nriter x nrt matrix of time points at which fT and FT are evaluated. Default value is [].
% qvec     - nriter x nrq matrix of quantiles at which qT is evaluated. Default value is [].
% tunnel   - Equals 0 (1) if stochastic tunneling is accounted for or not. Default value is 0.
%
% Outupt parameters
% ET     - nriter x 1 vector of expected values of waiting time T
% DT     - nriter x 1 vector of standard deviations of waiting time T
% fT     - length(tvec)xnriter matrix of values of the density function of T, at points in tvec.  
% FT     - length(tvec)xnriter matrix of values of the distribution function of T, at points in tvec.
% qT     - length(qvec)xnriter matrix of values of the quantile function of T, at points in qvec.

% Set default values of input parameters, if necessary
if nargin < 16
   tunnel = 0;
end
if nargin < 15
   qvec = [];
end
if nargin < 14
   tvec = [];
end

% Compute the number of iterations from the dimensionality of those input parameters that can be varied 
% between iterations. 
nriter = length(Deltat);
nriter = max(nriter,length(N));
nriter = max(nriter,length(L));
nriter = max(nriter,length(M));  
nriter = max(nriter,length(m));
nriter = max(nriter,length(mu));
nriter = max(nriter,length(gamma));
nriter = max(nriter,length(W));
nriter = max(nriter,length(K));
nriter = max(nriter,length(deltamax));
nriter = max(nriter,length(C));
nriter = max(nriter,length(tunnel));
ssize = size(s);
if (length(fitntype)==4) & (fitntype=='Mult')
   nriter = max(nriter,ssize(2));
elseif (length(fitntype)==5) & (fitntype=='Final')
   nriter = max(nriter,ssize(2));
elseif (length(fitntype)==5) & (fitntype=='Stair')
   nriter = max(nriter,ssize(1));
elseif (length(fitntype)==7) & (fitntype=='General')
   nriter = max(nriter,ssize(1));
end
%nriter

% Adjust the dimensionality of those input parameters that can be varied between iterations, so that they 
% are consistent with the number of iterations. If an input parameter has a single value, although the number of
% iterations is larger than one, this input parameter is kept constant.
if length(Deltat) < nriter
   Deltat = Deltat*ones(1,nriter);
end
if length(N) < nriter
   N = N*ones(1,nriter);
end
if length(L) < nriter
   L = L*ones(1,nriter);
end
if length(M) < nriter
   M = M*ones(1,nriter);
end
if length(m) < nriter
   m = m*ones(1,nriter);
end
if length(mu) < nriter
   mu = mu*ones(1,nriter);
end
if length(gamma) < nriter
   gamma = gamma*ones(1,nriter);
end
if length(W) < nriter
   W = W*ones(1,nriter);
end
if length(K) < nriter
   K = K*ones(1,nriter);
end
if length(deltamax) < nriter
   deltamax = deltamax*ones(1,nriter);
end
if length(C) < nriter
   C = C*ones(1,nriter);
end
if length(tunnel) < nriter
   tunnel = tunnel*ones(1,nriter);
end
if (length(fitntype)==4) & (fitntype=='Mult')
   if ssize(2) < nriter
      s = s*ones(1,nriter);
   end
elseif (length(fitntype)==5) & (fitntype=='Final')
   if ssize(2) < nriter
      s = s*ones(1,nriter);
   end
elseif (length(fitntype)==5) & (fitntype=='Stair')
   if ssize(1) < nriter
      s = ones(nriter,1)*s;
   end
elseif (length(fitntype)==7) & (fitntype=='General')
   if ssize(1) < nriter
      if ssize(2) < max(M)
         s(1,(ssize(2)+1):max(M)) = ones(1,max(M)-ssize(2));
      end
      s = ones(nriter,1)*s;
   end
end

% Compute waiting time statistics for all iterations
fT = [];
FT = [];
qT = [];
for i = 1:nriter
   if ((length(fitntype)==4) & (fitntype=='Mult')) | ((length(fitntype)==5) & (fitntype=='Final')) 
      [ET(i,1),DT(i,1),fTiter,FTiter,qTiter] = waitingtime_regseq_whaleevolution(Deltat(i),N(i),...
                                                  L(i),M(i),m(i),mu(i),gamma(i),W(i),fitntype,s(i),...
                                                  K(i),deltamax(i),C(i),tvec,qvec,tunnel(i)); 
   elseif (length(fitntype)==5) & (fitntype=='Stair') 
      [ET(i,1),DT(i,1),fTiter,FTiter,qTiter] = waitingtime_regseq_whaleevolution(Deltat(i),N(i),...
                                                  L(i),M(i),m(i),mu(i),gamma(i),W(i),fitntype,s(i,:),...
                                                  K(i),deltamax(i),C(i),tvec,qvec,tunnel(i)); 
   elseif (length(fitntype)==7) & (fitntype=='General')
      [ET(i,1),DT(i,1),fTiter,FTiter,qTiter] = waitingtime_regseq_whaleevolution(Deltat(i),N(i),...
                                                  L(i),M(i),m(i),mu(i),gamma(i),W(i),fitntype,...
                                                  s(i,1:(M(i)+1)),K(i),deltamax(i),C(i),tvec,qvec,tunnel(i));
   end
   fT = [fT fTiter];
   FT = [FT FTiter];
   qT = [qT qTiter];
end