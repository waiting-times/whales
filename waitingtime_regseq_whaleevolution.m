function [ET,DT,fT,FT,qT,kappa,Lambda,FT2,ET2] = waitingtime_regseq_whaleevolution(Deltat,N,L,M,m,mu,gamma,W,fitntype,s,K,deltamax,C,tvec,qvec,tunnel,nI,ETchoice)
% Compute various quantities of the waiting time until targeted binding sites have appeared due to
% mutations in the regulatory regions next to m out of M>=m possible genes. The order at which the targeted 
% binding sites at each gene appears is arbitrary, and for this reason the size of state space of the 
% the array component process C_t can be reduced from C^M to nrstates=M+1 if C=2 and (M+1)*(M+2)/2
% if C=3.  

% Input parameters:
% Deltat   - Generation time (typically in years)
% N        - Popoulation size
% L        - Length of regulatory region of each gene
% M        - Number of candidate genes to be part of the final target of m genes  
% m        - Number of genes that make up a target.
% mu       - mutation rate per generation and locus in each regulatory region
% gamma    - backward mutation probability
% W        - length of binding site
% fitntype - Specifies fitness function, either: 
%            'Mult' (multiplicative), 
%            'Final' (final target has different fitness)
%            'Stair' (intermediate and final types have different fitness)
%            'General' (general type of fitness funcrtion) 
% s        - Speciefies fitness parameters s(i) of an organism for which i=0,...,M local targets 
%            (binding sites) have been reached. The input format is either:
%            fitntype = 'Mult': s is positive scalar such that s(i) = s^i.
%            fitntype = 'Final': s is positive scalar such that s(i) = 1(i<m) + s*1(i>=m). 
%            fitntype = 'Stair': s is 1x2 vector [s1 s2] such that s(i) = 1(i=0) + s1*1(0< i < m) 
%                                + s2*1(i\ge m). 
%            fitntype = 'General': s is 1x(M+1)-vector [s(0) ... s(M)]. 
% K        - Number of possible binding sites per gene
% deltamax - Maximum number of mismatches that are allowed between a substring of length W of a 
%            regulatory sequence and the closest targeted binding site of that gene.
% C        - Number of intervals per gene for the distance that quantifies number of mismatches 
%            between a regulatory sequence and the closest targted binding site. These intervals are:
%            C = 2: {0,...,deltamax},{deltamax,...,W}
%            C = 3: {0,...,deltamax},{deltamax+1},{deltamax+2,...,W}
% tvec     - Vector of time points at which FT and fT are evaluated. Default value is [].
% qvec     - Vector of quanties at which qT is evaluated. Default value is [].
% tunnel   - Equals 1 (0) if stochastic tunneling is accounted for or not. Default value is 0.
% nI       - A positive integer that corresponds to the number of intervals used for computing FT2 
%            through numerical integration of the density function. The default value is nI=[].
%            Then FT2 is not computed.
% ETchoice - Equals 1 (0) if an alternative explicit formula for ET is computed or not. This formul
%            is only applicable however for C=2 and fitntype='Mult'. Default value is 0.

% Output paramters
% ET     - Expected value of waiting time T
% DT     - Standard deviation of waiting time T
% fT     - length(tvec)x1 vector with density function of T evaluated at points in tvec.  
% FT     - length(tvec)x1 vector with distribution function of T evaluated at points in tvec.
% qT     - length(qvec)x1 vector with quantile function of T evaluated at points in qvec.
% kappa  - nrstates x 1 vector with initial distribution of C_t
% Lambda - nrstates x nrstates intensity matrix of C_t 
% FT2    - length(tvec)x1 vector with an alternative way of computing the distribution function of
%          T at the time poinst stored in tvec. This computation is based on numberical integration 
%          of the density function, and method that is preferable when FT is very small (since it 
%          avoids cancellation when two numbers, almost equally large, are subtracted, as in the 
%          formula for FT).
% ET2    - Alternative expression for ET, based on an explicit formula. It is only computed if
%          ETchoice=1, C=2 and fitntype='Mult'.

% Set default values of input parameters, if necessary
if nargin < 18
   ETchoice = 0;
end
if nargin < 17
   nI = [];
end
if nargin < 16
   tunnel = 0;
end
if nargin < 15
   qvec = [];
end
if nargin < 14
   tvec = [];
end

% When C=3, use the fact that fitntype='Final' is a special case of fitntype='General'.
if (C == 3) & ((length(fitntype)==5) & (fitntype=='Final')) 
   fitntype = 'General';
   sfinal = s;
   s = ones(1,M+1);
   s((m+1):(M+1)) = sfinal;
end

% Use the fact tht fitntype='Stair' is a special case of fitntype='General'.
if ((length(fitntype)==5) & (fitntype=='Stair')) 
   fitntype = 'General';
   sintermediate = s(1);
   sfinal = s(2);
   s = ones(1,M+1);
   s(2:m) = sintermediate;
   s((m+1):(M+1)) = sfinal;
end

% Find the stationary distribution and transition rates of each gene
if C == 2
    if (length(fitntype)==5) & (fitntype=='Final')
      [ETone,DTone,fTone,FTone,qTone,kappaone,Lambdaone] ...
          = waitingtime_regseq_2(N,L,1,1,gamma,W,[1 1],K,[0 deltamax+1],deltamax,'arbitrary',[],[],[],0);
    elseif (length(fitntype)==4) & (fitntype == 'Mult')
      [ETone,DTone,fTone,FTone,qTone,kappaone,Lambdaone] ...
          = waitingtime_regseq_2(N,L,1,1,gamma,W,[1 s],K,[0 deltamax+1],deltamax,'arbitrary',[],[],[],0);
    elseif (length(fitntype)==7) & (fitntype == 'General')
       for i = 1:M
          [ETone,DTone,fTone,FTone,qTone,kappaone,Lambdaone(:,:,i)] ...
          = waitingtime_regseq_2(N,L,1,1,gamma,W,[1 s(i+1)/s(i)],K,[0 deltamax+1],deltamax,'arbitrary',[],[],[],0);
       end  
    end
elseif C == 3
   if (length(fitntype)==4) & (fitntype == 'Mult')
      [ETone,DTone,fTone,FTone,qTone,kappaone,Lambdaone] ...
         = waitingtime_regseq_2(N,L,1,1,gamma,W,[1 s],K,[0 deltamax+1 deltamax+2],deltamax,'arbitrary',[],[],[],0);
   elseif (length(fitntype)==7) & (fitntype == 'General')
      for i = 1:M
         [ETone,DTone,fTone,FTone,qTone,kappaone,Lambdaone(:,:,i)] ...
            = waitingtime_regseq_2(N,L,1,1,gamma,W,[1 s(i+1)/s(i)],K,[0 deltamax+1 deltamax+2],deltamax,'arbitrary',[],[],[],0);
      end
   end
end
%ETone
%kappaone
%Lambdaone

% Compute number of states of the array component process C_t
if C == 2
   nrstates = M+1;
elseif C == 3
   nrstates = (M+1)*(M+2)/2;
end

% If C=3, create list of codings of all states (j,k), where 0 <= j <= M refers to the number of 
% genes one distance class away from the local target, whereas 0 <= k <= M refers to the number of
% genes that have reached the target, so that j+k <= M. i=order(j+1,k+1) gives the order number of 
% state (j,k) whereas state(i,:) = [j k] gives the state that corresponds to order number i. Finally
% compute (for C=2 and C=2) the number of genes g(i) that have reached the local target, for each 
% state i.
if C == 2
   g = 0:M;
elseif C == 3
   order = zeros(M+1,M+1);
   state = [];
   i = 0;
   for k = 0:M
      for j = 0:(M-k)
         i = i+1;
         order(j+1,k+1)=i;
         state = [state; j k];
         g(i) = k;
      end
   end
end
%g

% Compute kappa, the initial distribution of the array component process C_t at time t=0.
if C == 2
   kappa = binopdf(0:M,M,kappaone(2));
elseif C == 3
   for i = 1:nrstates
      j = state(i,1);
      k = state(i,2);
      kappa(i) = nchoosek(M,k)*nchoosek(M-k,j) ...
               * kappaone(1)^(M-j-k)*kappaone(2)^j*kappaone(3)^k;
   end
end

% Compute Lambda, a matrix with transition rates between all paira of states for the array component
% process C_t
Lambda = zeros(nrstates,nrstates);
if C == 2
   if (length(fitntype)==7) & (fitntype == 'General')
      for i = 1:M
         Lambda(i,i+1) = (M-i+1)*Lambdaone(1,2,i);
         Lambda(i+1,i) = i*Lambdaone(2,1,i);
      end
   else  
      for i = 1:M
         Lambda(i,i+1) = (M-i+1)*Lambdaone(1,2);
         Lambda(i+1,i) = i*Lambdaone(2,1);
      end
      if (length(fitntype)==5) & (fitntype == 'Final')
         if s ~= 1
            Lambda(m,m+1) = Lambda(m,m+1)*N*(1-s^(-1))/(1-s^(-N)); % Adjust rates for when fitness type is 'Final'.
            if m < M
               Lambda(m+1,m) = Lambda(m+1,m)*N*(1-s)/(1-s^N);
            end
         end
      end
   end
   Lambda(1,1) = -Lambda(1,2);
   for i = 2:M
      Lambda(i,i) = - (Lambda(i,i-1)+Lambda(i,i+1));
   end
   Lambda(M+1,M+1) = - Lambda(M+1,M);
elseif C == 3
   if (length(fitntype)==7) & (fitntype == 'General')
      for i = 1:nrstates
         j = state(i,1);
         k = state(i,2);
         if j > 0
            idown = order((j-1)+1,k+1);
            Lambda(i,idown) = j*Lambdaone(2,1,k+1);
            iup = order((j-1)+1,(k+1)+1);
            Lambda(i,iup) = j*Lambdaone(2,3,k+1);
         end
         if (j+k) < M
            iup = order((j+1)+1,k+1);
            Lambda(i,iup) = (M-k-j)*Lambdaone(1,2,k+1);
         end
         if k > 0
            idown = order((j+1)+1,(k-1)+1);
            Lambda(i,idown) = k*Lambdaone(3,2,k);
         end
      end
   elseif (length(fitntype)==4) & (fitntype == 'Mult')
      for i = 1:nrstates
         j = state(i,1);
         k = state(i,2);
         if j > 0
            idown = order((j-1)+1,k+1);
            Lambda(i,idown) = j*Lambdaone(2,1);
            iup = order((j-1)+1,(k+1)+1);
            Lambda(i,iup) = j*Lambdaone(2,3);
         end
         if (j+k) < M
            iup = order((j+1)+1,k+1);
            Lambda(i,iup) = (M-k-j)*Lambdaone(1,2);
         end
         if k > 0
            idown = order((j+1)+1,(k-1)+1);
            Lambda(i,idown) = k*Lambdaone(3,2);
         end
      end
   end
   for i = 1:nrstates
      if i > 1  
         Lambda(i,i) = - sum(Lambda(i,1:(i-1)));
      end
      if i < nrstates
         Lambda(i,i) = Lambda(i,i) - sum(Lambda(i,(i+1):nrstates));
      end
   end
end
Lambda = Lambda*mu/Deltat;

% Compute states that correspond to absorbing and non-absorbing states, then define the submatrices
% of Lambda corresonding to transitions among non-absorving states and from non-absorbing to absorbing
% states respectively.
Inonabs = find(g<m);
Iabs = find(g>=m);
Lambda0 = Lambda(Inonabs,Inonabs);
lambda1 = Lambda(Inonabs,Iabs)*ones(length(Iabs),1);
kappa0 = kappa(Inonabs);
kappanonabs = sum(kappa0);

% The wating time is defined as the time when the array component proccess C(t) reaches an 
% absorbing state. Compute the expected value, standard deviation, density function, distribution 
% function and quantiles of the waiting time distribution, assuming that C(t) is a continuous time 
% Markov process, so T has a phase-type distribution.
onevec = ones(length(Inonabs),1);
ET = - kappa0*inv(Lambda0)*onevec;
DT = sqrt(2*kappa0*inv(Lambda0)^2*onevec-ET^2);
if length(tvec)==0
      fT = [];
      FT = [];
else
   for i=1:length(tvec)
      fT(i,1) = kappa0*expm(Lambda0*tvec(i))*lambda1;
      FT(i,1) = 1 - kappa0*expm(Lambda0*tvec(i))*onevec;
   end
end
if length(qvec)==0
   qT = [];
else
   if length(Inonabs)==1 % Waiting time has an exponential distribution if there is only one non-absorbing state.
      qT = (log(min([ones(1,length(qvec));(1-qvec)/kappa0])))'/Lambda0;
   else
      for i=1:length(qvec)
         if qvec(i)<=1-kappanonabs
            qT(i,1) = 0;
         else 
            qTint = [0,ET+50*DT];
            qTmid = (qTint(1)+qTint(2))/2;
            iter = 0;
            while (((qTint(2)-qTint(1))/qTmid) > 10^(-6)) & (iter<100)
               qTmid = (qTint(1)+qTint(2))/2;
               F = 1 - kappa0*expm(Lambda0*qTmid)*onevec;
               if F > qvec(i)
                  qTint(2)=qTmid;
               else 
                  qTint(1)=qTmid;
               end
%               qTint
               iter = iter+1;
            end
            if iter < 100
               qT(i,1) = (qTint(1)+qTint(2))/2;
            else 
               qT(i,1) = Inf;
            end
         end
      end
   end
end

% Alternative computation FT2 of the distribution function, based on numerical integration of
% the density function. This numerical integration uses Simpson's formula over G intervals.
if length(nI) == 0
   FT2 = [];
else
   nG = 2*nI+1;  % Number of grid poinst when Simpson's formula is applied to nI intervals.
   for i=1:length(tvec)
      if tvec(i) == 0
         FT2(i,1) = kappa(nrstates); % Probability that T equals 0.
      else
         tvecint = 0:(tvec(i)/(2*nI)):tvec(i);
         for j = 1:nG
            fTint(j) = kappa0*expm(Lambda0*tvecint(j))*lambda1;
         end
         FT2(i,1) = (sum(fTint([1 nG])) + 2*sum(fTint(3:2:(nG-2))) + 4*sum(fTint(2:2:(nG-1))))/(6*nI);
         FT2(i,1) = FT2(i,1) * tvec(i);
         FT2(i,1) = FT2(i,1) + kappa(nrstates); % Add probability that T equals 0.
      end
   end
end

% Alternative explicit expression ET2 for the expected waiting time.
if (ETchoice==1) & (C==2) & (length(fitntype)==4) & (fitntype=='Mult')
   forwrate = Lambdaone(1,2)*mu/Deltat; % Forward rate of one gene.
   r = Lambdaone(2,1)/Lambdaone(1,2);  % Ratio between backward and forward rates for one gene.
   ET2 = 0;
   for i = 1:m
      theta = 0;
      for h = (i-1):(m-1)
         for k = 0:h
            theta = theta + nchoosek(M-1,k)*r^(h-k)/((M-k)*nchoosek(M-1,h));
         end
      end
      ET2 = ET2 + kappa(i)*theta; 
   end 
   ET2 = ET2/forwrate;
else
   ET2 = [];     
end

% Final adjustment of output parameters
kappa = kappa';

