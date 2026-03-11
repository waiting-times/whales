function [ET,DT,fT,FT,qT,kappa,Lambda,Inonabs,Lambda0,g,cmatrix,EH,VH,M,Meta] = waitingtime_regseq_2(N,L,m,mu,gamma,W,s,K,delta,deltamax,TA,tvec,qvec,b,tunnel)
% Consider a haploid population of size N that evolves in continuous time according to a Moran model
% with genetic drift, mutations, and selection. To each individual we associate m distinct
% regulatory sequences of length L nucleotides, written in a four letter alphabet 1=A, 2=C, 3=G, 4=T.
% These m sequences are located in the promoter regions of m distinc genes. It is assumed that the 
% mutation probability mu per generation, nucleotide, and individual is small (mu << 1/(NL)). 
% Genetic drift will then dominate over mutations, so that most individuals at any time point 
% have almost the same mxL array of regulatory sequences within the m promoter regions. This 
% consensus array of dimension mxL will vary over time due to mutations, drift and selection, and T is 
% the waiting time until m pre-specified targets have appeared in each one of the m regulatory regions.
% The target of regulatory sequence j (j=1,...,m) is defined in such a way that some substring of the 
% regulatory sequence of length W (2\le W < L) should be within distance deltamax to one of K(j)
% pre-specified bindings sites b_{j1},...,b_{jK(j)} of length W. It is further assumed that whenever
% a nucleotide mutates, it is equally likely (=1/3) that any of the other three nucleotides appears
% after the mutation.  
%
% We approximate the distribution of T by a phase-type distribution, as described in Hoessjer, Bechly, and 
% Gauger (2018). More specifically, all regulatory arrays are defined into a finite number of components,
% and the process C_t=(C_{t1},...,C_{tm}; t>=0) monitors which components the arrays belongs to at time t. 
% It is assumed that C_t is a continuous time Markov process, and T is the time when C_t first hits an 
% absorbing state.
%
% Compared to version 1 (waitingtime_regseq), this version 2 (waitingtime_regseq_2) has an improved/corrected
% formula for the way in which stochastic tunneling affect the intensity matrix of the Markov process C_t. 
%
% Input parameters:
% N        - population size
% L        - length of each regulatory sequence
% m        - number of promoter regions
% mu       - mutation probability
% gamma    - back mutation probability
% W        - length of binding sites 
% s        - 1x(m+1) vector of selection coefficients, where s(j) is the relative fitness of a regulatory
%            array for which j-1 targets have appeared, compared to an array for which none of the targets
%            has appeared. Without loss of generality, one lets s(1)=1.
% K        - 1xm vector of number of possible binding site targets at the different regulatory sequences.
%            If K is a scalar and m>1, then K(j)=K for j=1,...,m.
% delta    - 1xC vector of left end points of C intervals (c=1,...,C) that for each regulatory region
%            specifies how far away it is from the closest binding site of that region. This is quantified
%            in terms of the number of nucleotides that differ from the closes binding site, and interval
%            c consists of distances delta(C-c+1),...,delta(C-c+2)-1, where 0 = delta(1) < delta(2) < ...
%            < delta(C)<=W. 
% deltamax - maximal distance of a regulatory sequence, to the closest binding site of that 
%            promoter region, in order to have reached the target of this region.
% TA       - string that specifies type of target appearanece. It equals 'fixed' or 'arbitrary',
%            depending on whether the targets have to appear in a pre-specified order 1,...,m, or
%            could appear in any order.
% tvec     - vector of waiting times for which the distribution function FT and density function fT
%            of T are evaluated. If tvec = [] (default choice), then FT and fT are not computed. 
% qvec     - vector of quantiles of the distribution of T that are evaluated (=qT). If qvec = [] (default choice), 
%            then qT is not evaluated.
% b        - nrbxW-matrix that specifies the targeted binding sites. 
%            If nrb=sum(K)-m, then no overdispersion is assumed for the number of binding site hits along
%            the regulatory sequences of the m genes. We may then assume, without loss of generality, 
%            that the first binding site of region j is [1...1]. The other K(j)-1 binding
%            sites of region j are specified as rows K(1)+...+K(j-1)-(j-1)+1,...,K(1)+...+K(j)-j of b.
%            In particular, b is omitted (b=[]) if K(1)=...=K(m)=1.
%            If b is omitted (b=[]), even though some gene j has K(j)>1, an approximate algorithm is 
%            used for calculating the probability vector M(j,:) and the intensity matrices M1(:,:,j)
%            and M2(:,:,j). This approximate algorithm assumes that the b(j) binding sites of 
%            regulatory region j are so widely dispersed that subwords of the regulatory sequence has 
%            a closest mismatch only to one binding site of that gene, whenever the elements of kappa 
%            and M are calculated. If K(j) is too large, there will inevitably be overlaps beween the
%            "mismatch neighborhoods" of the K(j) binding sites. 
%            If nrb=sum(K), then overdispersion is acounted for in the number of binding site hits along 
%            the regulatory sequences of the m genes. The K(j) binding sites of region j are specified
%            as rows K(1)+...+K(j-1)+1,...K(1)+...+K(j) of b.
% tunnel   - Indicator of whether stochastic tunneling is included (tunnel=1) or not (tunnel = 0).
%            Default choice of tunnel is 1. 
%
%
% Output parameters:
% ET      - expected value of T
% DT      - standard deviation of T
% fT      - length(t)x1 vector of values of the density function of T, at the time points specified in t.
% FT      - length(t)x1 vector of values of the distribution function of T, at the time points specified in t.
% qT      - length(t)x1 vector of quantiles of T, at quantiles specified by q.
% kappa   - (C^m)x1 vector, corresponding to the initial distribution of the array component process C_t. 
% Lambda  - (C^m)x(C^m)-matrix, corresponding to the transition matrix of the array component process C_t (before any states have been identified as absorbing).
% Inonabs - 1xlength(Inonabs) vector, containing the indeces (in {1,...,C^m}) of all nonabsorbing states. 
% Lambda0 - length(Inonabs)xlength(Inonabs) matrix - a submatrix of Lambda with transition rates between all non-absorbing states.
% g       - (C^m)x1 vector containing the types of array components
% cmatrix - (C^m)xm matrix containing all elements of c=[c(1),...,c(m)] of the array component process C_t, listed in lexicographc order
% EH      - mxC matrix, where EH(j,c) gives the expected number of words of length W in region j that
%           belong to distance class c to the closest binding site of that region.
% VH      - mxC matrix, where VH(j,c) gives the variance of the number of words of length W in region j 
%           that belong to distance class c to the closest binding site of that region.
% M       - mxC matrix, where M(j,c) gives the probability the a randomly chosen word of length W in region j
%           belongs to distance class c.
% Meta    - mxCx(W-1) array, where Meta(j,c,eta) gives the probability that two randomly chosen words 
%           of length W both belong to distance class c, given that the last W-eta letters of word 1
%           agree with the last W-eta letters of word 2. This output variable is only computed when 
%           overdispersion is accounted for, otherwise Meta = []. 

% Initializations
C = length(delta); % Number of possible distance intervals between a regulatory seq and target.
delta = [delta W+1]; % Augment delta so that right end point (=W) of interval C is well defined.
find(delta==(deltamax+1));
Cpr = C+1-find(delta==(deltamax+1)); % Number of distance intervals for which target has not been reached.
intnr = zeros(1,W+1); % Intnr is a 1x(W+1) vector, where intnr(i) specifies which interval contains i-1.
                      % This is needed in order to map a distance to the closest binding site, to a binding site class.
for i=1:W+1
   intnr(i) = C+1-sum(delta<=(i-1));
end
%intnr
if length(K)==1 % If K is a scalar, then K(j)==K.
   K = K*ones(1,m);
end
L0 = L-W+1; % Number of possible starting points of a binding site of length W along a regulatory region of length L.
if nargin<12
   tvec = []; % Defalut choice (not evaluating FT and fT)
end
if nargin<13
   qvec = []; % Default choice (not evaluating qT)
end
if nargin<14
   overdisp = 0; % No overdispersion when b is not specified
elseif nargin >= 14
   [nrb nrletters] = size(b);
   if nrb == sum(K)
      overdisp = 1; % Overdispersion when all K(j) binding sites are specified at all m genes.
   elseif nrb == sum(K)-m
      overdisp = 0; % No overdispersion when only K(j)-1 binding sites per gene are specified.
   elseif (nrb == 0) | (sum(K)>1)
      overdisp = 0; % No overdispersion when some K(j)>1 and the approximate algorithm is used.
   end
end
if nargin<15
   tunnel = 1; % Default choice of including stochastic tunneling. 
end

% Compute vector g = [g(1),...,g(C^m)], where g(i) gives the number of reached targets
% for an array component that is indexed by a vector c=[c(1),...,c(m)], with 1<=c(j)<=C.
% It is assumed that these vectors are ordered lexicographically, so that i is the order number of
% vector c. A target is reached within promoter (regulatory) region j if c(j)>Cpr.
cmatrix = [];
if (length(TA)==5) & (TA=='fixed')
   for i = 1:C^m
      c = invlexorder(i,C,m);
      cmatrix = [cmatrix; c];
      if min(c)<=Cpr
         g(i) = min(find(c<=Cpr))-1;
      else
         g(i) = m;
      end
   end
elseif (length(TA)==9) & (TA=='arbitrary')
   for i = 1:C^m
      c = invlexorder(i,C,m);
      cmatrix = [cmatrix; c];
      g(i) = sum(c>Cpr);
   end
end
%g
%s(g+1)

% Compute set of non-absorbing and absorbing states. They are both subsets of {1,...,C^m}.
Inonabs = find(g<m);
Iabs = find(g==m);

% Compute the mxC matrix M=(M(j,c)). Each row of M is a probability distribution, where M(j,c) is the 
% probability that a randomly chosen word of length W has a distance to the closest binding site of 
% promoter region j that belongs to interval number c (c=1,...,C).
M = zeros(m,C);
for j=1:m
   if (K(j)==1) | ((K(j)>1) & (length(b)==0)) % K(j)>1 corresponds to an approximate algorithm
      for c=1:C
         for a=delta(C-c+1):(delta(C-c+2)-1)
            M(j,c) = M(j,c) + nchoosek(W,a)*3^a/4^W;
         end
      end
      if K(j) > 1 % The approximate algorithm for K(j)>1 (used when b=[])
         M(j,2:C) = K(j)*M(j,2:C);
         M(j,1) = 1 - sum(M(j,2:C));
         if M(j,1) < 0 % When K(j) is too large for K(j) binding sites to be widely dispersed, an error message (NaN) is output.
            M(j,1) = NaN;
         end
      end          
   elseif K(j)==2
      if overdisp == 0
         b1 = ones(1,W);
         b2 = b(sum(K(1:j))-j,:);
      elseif overdisp == 1
         b1 = b(sum(K(1:j))-1,:);
         b2 = b(sum(K(1:j)),:);
      end
      Wj = sum(b1~=b2);
      for v1=0:Wj
         for v2=0:Wj-v1
            for c=1:C
               v3lower = max(0,delta(C-c+1)-Wj+max(v1,v2));
               v3upper = min(delta(C-c+2)-1-Wj+max(v1,v2),W-Wj);
               if v3upper >=0
                  for v3=v3lower:v3upper
                     M(j,c) = M(j,c) + prod(1:Wj)/(prod(1:v1)*prod(1:v2)*prod(1:Wj-v1-v2))*nchoosek(W-Wj,v3)*2^(Wj-v1-v2)*3^v3; 
                  end
               end
            end
         end
      end
      M = M/4^W;
   else
%      target = ones(1,W); % This is used for testing the general algorithm when K(j)=1. 
      if overdisp == 0
         target = [ones(1,W); b(sum(K(1:j-1))-(j-1)+1:sum(K(1:j))-j,:)];
      elseif overdisp == 1
         target = b(sum(K(1:j-1))+1:sum(K(1:j)),:);
      end
      for i=1:4^W  % Go through all 4^W words of length W and compute distance to closest binding site.
         word = invlexorder(i,4,W);
         distance = min(((ones(K(j),1)*word)~=target)*ones(W,1));
         M(j,intnr(distance+1)) = M(j,intnr(distance+1))+1;
      end
      M(j,:) = M(j,:)/4^W;
   end
end

% Compute the mxCx(W-1) array Meta=(Meta(j,c,eta)), where Meta(j,c,eta) is the probability that a
% randomly chosen pair of words of length W both belong to distance class c, given that the first 
% W-eta letters of the first word agree with the last W-eta letters of second word.
if overdisp == 0
   Meta = [];
elseif overdisp == 1
   for eta=1:W-1
      Meta(:,:,eta) = zeros(m,C);
      for i=1:4^(W-eta)  % Go through all 4^(W-eta) subwords of length W-eta that are identical for the two words.
         subword1 = invlexorder(i,4,W-eta);
         for j=1:m
            target = b(sum(K(1:j-1))+1:sum(K(1:j)),:);
            sum1 = zeros(1,C);
            sum2 = zeros(1,C);
            for k=1:4^eta
               subword2 = invlexorder(k,4,eta);
               word1 = [subword1 subword2];
               distance1 = min(((ones(K(j),1)*word1)~=target)*ones(W,1));
               sum1(intnr(distance1+1)) = sum1(intnr(distance1+1))+1;
               word2 = [subword2 subword1];
               distance2 = min(((ones(K(j),1)*word2)~=target)*ones(W,1));
               sum2(intnr(distance2+1)) = sum2(intnr(distance2+1))+1;
            end
            for c=1:C
               Meta(j,c,eta) = Meta(j,c,eta)+sum1(c)*sum2(c); 
            end        
         end
      end
      Meta(:,:,eta) = Meta(:,:,eta)/4^(W+eta);
   end
end

% Compute the mxC matrices EH=(EH(j,c)) and VH=(VH(i,j), where EH(j,c) is the expected number of words
% of length W along the regulatory sequence of promoter region j, whose distance to the closest 
% binding site of that region belongs to distance class c (c=1,...,C). This is also referred to as 
% the expected number of words that hit distance class c for region j. Similarly, VH(j,c) is the 
% variance of this hit variable. 
EH = L0*M;
if overdisp == 0 % A Poisson distributed number of hits is assumed for all hit variables.
   VH = EH;
elseif overdisp == 1 % A Poisson or negative binomial distrribution is assumed for the hit variables when
                     % overdispersion is allowed for, depending on whether there is overdispersion or not.
   VH = L0*(M-M.^2);
   for j=1:m
      for c=1:C
         for eta=1:(W-1)
            VH(j,c) = VH(j,c) + 2*(L0-eta)*(Meta(j,c,eta)-M(j,c)^2);
         end
      end
   end
end

% Compute the CxCxm array M1=(M1(c,d,j)). For each j, M1(.,.,j) is a matrix, containing the single 
% mutation probability flows M1(c,d,j) between all pairs Bjc and Bjd of sets of words of length W. 
% Here Bjc and Bjd correspond to all words of length W  that have a distance to the closest targeted 
% binding site of promoter region j that belongs to intervals number c and d respectively. A 
% transition from Bjc to Bjd occurs when a single letter of a word in Bjc is changed, such that the 
% new word belongs to Bjd. 
for j=1:m
   M1(:,:,j) = zeros(C,C);
   if (K(j)==1) | ((K(j)>1) & (length(b)==0)) % K(j)>1 corresponds to an approximate algorithm
      for c=1:C
          d = c+1;
          if d<=C
             M1(c,d,j) = K(j)*nchoosek(W,delta(C-c+1))*3^delta(C-c+1)*delta(C-c+1)/(3*4^W);
          end
          d = c-1;
          if d>=1
             M1(c,d,j) = K(j)*nchoosek(W,delta(C-c+2)-1)*3^delta(C-c+2)*(W-delta(C-c+2)+1)/(3*4^W);
          end
      end
   elseif K(j)==2
      if overdisp == 0
         b1 = ones(1,W);
         b2 = b(sum(K(1:j))-j,:);
      elseif overdisp == 1
         b1 = b(sum(K(1:j))-1,:);
         b2 = b(sum(K(1:j)),:);         
      end
      Wj = sum(b1~=b2);
      for v1=0:Wj
         for v2=0:Wj-v1
            for c=1:C-1
               v3 = delta(C-c+1)-Wj+max(v1,v2);
               if (v3>=0) & (v3<=W-Wj)
                  M1term = prod(1:Wj)/(prod(1:v1)*prod(1:v2)*prod(1:Wj-v1-v2))*nchoosek(W-Wj,v3)*2^(Wj-v1-v2)*3^v3; 
                  if v1>v2
                     M1term = M1term*(Wj+v3-v1);
                  elseif v1==v2
                     M1term = M1term*(2*Wj+v3-v1-v2);
                  elseif v1<v2
                     M1term = M1term*(Wj+v3-v2);
                  end
                  M1(c,c+1,j) = M1(c,c+1,j) + M1term;
               end
            end
         end
      end
      for c=2:C
         M1(c,c-1,j) = M1(c-1,c,j);
      end
      M1(:,:,j) = M1(:,:,j)/(3*4^W);
   else
%      target = ones(1,W); % This is used for testing the general algorithm when K(j)=1.
      if overdisp == 0
         target = [ones(1,W); b(sum(K(1:j-1))-(j-1)+1:sum(K(1:j))-j,:)];
      elseif overdisp == 1
         target = b(sum(K(1:j-1))+1:sum(K(1:j)),:);
      end
      for i=1:4^W % Go through all 4^W words (=word) of length W and compute distance to closest binding site. 
         word = invlexorder(i,4,W);
         distance = min(((ones(K(j),1)*word)~=target)*ones(W,1));
         c = intnr(distance+1); % Interval number that the distance  of word belongs to.
         for w=1:W % For all 3*W neighbouring words (=neighword) of word, compute distance to closest binding site.
            for a=1:4
               if word(w)~=a
                  neighword = word; 
                  neighword(w) = a;
                  neighdistance = min(((ones(K(j),1)*neighword)~=target)*ones(W,1));
                  d = intnr(neighdistance+1); % Interval number that the distance of the neighbouring word belongs to.
                  if c~=d
                     M1(c,d,j) =  M1(c,d,j)+1;
                  end
               end
            end
         end
      end
      M1(:,:,j) = M1(:,:,j)/(3*4^W); % Convert number of neighbouring pairs of words that belong
                                     % to different distance classes, to probability flow unit.
   end
end
%M1

% Compute the CxCxm array M2=(M2(c,d,j)). For each j, M2(.,.,j) is a matrix, containing the double 
% mutation probability flows M2(c,d,j) between all pairs Bjc and Bjd of sets of words of length W. 
% Here Bjc and Bjd correspond to all words of length W  that have a distance to the closest targeted 
% binding site of promoter region j that belongs to intervals number c and d respectively. A 
% transition from Bjc to Bjd occurs when a two letters of a word in Bjc are changed, such that the 
% new word belongs to Bjd, AND Bjd is a least two distance classes away from Bjc (|d-c|>=2). 
for j=1:m
   M2(:,:,j) = zeros(C,C);
   if (K(j)==1) | ((K(j)>1) & length(b)==0)  % K(j)>1 corresponds to an approximate algorithm
      for c=1:C
          d = c+2;
          if (d<=C) & (delta(C-c)==(delta(C-c+1)-1))
             M2(c,d,j) = K(j)*nchoosek(W,delta(C-c+1))*3^delta(C-c+1)*nchoosek(delta(C-c+1),2)/(9*4^W);
          end
          d = c-2;
          if (d>=1) & (delta(C-c+3)==delta(C-c+2)+1)
             M2(c,d,j) = K(j)*nchoosek(W,delta(C-c+2)-1)*3^(delta(C-c+2)+1)*nchoosek(W-delta(C-c+2)+1,2)/(9*4^W);
          end
      end
   elseif K(j)==2
      if overdisp == 0
         b1 = ones(1,W);
         b2 = b(sum(K(1:j))-j,:);
      elseif overdisp == 1
         b1 = b(sum(K(1:j))-1,:);
         b2 = b(sum(K(1:j)),:);         
      end
      Wj = sum(b1~=b2);
      for v1=0:Wj
         for v2=0:Wj-v1
            for c=1:C-2
               if delta(C-c)==(delta(C-c+1)-1)
                  v3 = delta(C-c+1)-Wj+max(v1,v2);
                  if (v3>=0) & (v3<=W-Wj)
                     M2term = prod(1:Wj)/(prod(1:v1)*prod(1:v2)*prod(1:Wj-v1-v2))*nchoosek(W-Wj,v3)*2^(Wj-v1-v2)*3^v3; 
                     if v1>v2
                        M2term = M2term*(Wj+v3-v1)*(Wj+v3-v1-1)/2;
                     elseif v1==v2
                        M2term = M2term*((Wj+v3)*(Wj+v3-1)/2-v1*v2+v3*(Wj-v1-v2)+(Wj-v1-v2)*(Wj-v1-v2-1)/2);
                     elseif v1<v2
                        M2term = M2term*(Wj+v3-v2)*(Wj+v3-v2-1)/2;
                     end
                     M2(c,c+2,j) = M2(c,c+2,j) + M2term;
                  end
               end
            end
         end
      end
      for c=3:C
         M2(c,c-2,j) = M2(c-2,c,j);
      end
      M2(:,:,j) = M2(:,:,j)/(9*4^W); % Convert number of neighbouring pairs of words that belong
                                     % to different distance classes, to probability flow unit.
   else
 %     target = ones(1,W); % This is used for testing the general algorithm when K(j)=1.
      if overdisp == 0
         target = [ones(1,W); b(sum(K(1:j-1))-(j-1)+1:sum(K(1:j))-j,:)];
      elseif overdisp == 1
         target = b(sum(K(1:j-1))+1:sum(K(1:j)),:);
      end
      for i=1:4^W % Go through all 4^W words (=word) of length W and compute distance to closest binding site. 
         word = invlexorder(i,4,W);
         distance = min(((ones(K(j),1)*word)~=target)*ones(W,1));
         c = intnr(distance+1); % Interval number that the distance  of word belongs to.
         for w1=1:W % For all 9*W*(W-1)/2 strings two letters from word (=neigh2word), compute distance to closest binding site.
            for w2=(w1+1):W           
               for a1=setdiff(1:4,word(w1))
                  for a2 = setdiff(1:4,word(w2))
                      neigh2word = word; 
                      neigh2word(w1) = a1;
                      neigh2word(w2) = a2;
                      neigh2distance = min(((ones(K(j),1)*neigh2word)~=target)*ones(W,1));
                      d = intnr(neigh2distance+1); % Interval number that the distance of neigh2word belongs to.
                      if abs(d-c)>= 2
                         M2(c,d,j) =  M2(c,d,j)+1;
                      end
                  end
               end
            end
         end
      end
      M2(:,:,j) = M2(:,:,j)/(9*4^W); % Convert number of neighbouring pairs of words that belong
                                     % to different distance classes, to probability flow unit.
   end
end
%M2

% Compute the initial distribution kappa = [kappa(1),...,kappa(C^m)]' of the array components process 
% C(t). kappa0 is the corresponding vector of initial probabilities for the non-absorbing states.
kappa = zeros(C^m,1);
for i = 1:C^m
   c = invlexorder(i,C,m);
   kappa(i) = 1;
   for j=1:m
      if overdisp == 0
         kappa(i) = kappa(i)*(1-exp(-EH(j,c(j))));
         if c(j)<C
            kappa(i)=kappa(i)*exp(-sum(EH(j,(c(j)+1):C)));
         end
      elseif overdisp == 1
         if VH(j,c(j)) <= EH(j,c(j))
            kappa(i) = kappa(i)*(1-exp(-EH(j,c(j))));
         else 
            a = EH(j,c(j))^2/(VH(j,c(j))-EH(j,c(j)));
            kappa(i) = kappa(i)*(1-(1+EH(j,c(j))/a)^(-a));
         end
         if c(j)<C
            for d=(c(j)+1):C
               if VH(j,d) <= EH(j,d)
                  kappa(i) = kappa(i)*exp(-EH(j,d));
               else
                  a = EH(j,d)^2/(VH(j,d)-EH(j,d));
                  kappa(i) = kappa(i)*(1+EH(j,d)/a)^(-a); 
               end
            end
         end
      end
   end
end
kappa=kappa/sum(kappa); % Normalize to be a probability distribution. 
kappa0 = kappa(Inonabs');
kappanonabs = sum(kappa0);

% Compute transition matrix Lambda=(Lambda(i,k)) of array component process C(t).
Lambda = zeros(C^m,C^m);
for i=1:C^m
   c = invlexorder(i,C,m);
   for j=1:m
      for step = [-1 1]
         d = c;
         d(j) = c(j)+step;
         if (max(d)<=C) & (min(d)>=1)
            k = lexorder(d,C,m);
            Lambda(i,k) = N*mu*M1(c(j),d(j),j)/M(j,c(j)); % This is transition rate for a neutral model                                                            % without backmutations when d(j)=c(j)+1.
            if overdisp == 0
               if step == 1
                  Lambda(i,k) = Lambda(i,k)*EH(j,c(j));
               elseif step == -1
                  Lambda(i,k) = Lambda(i,k)*EH(j,c(j))*exp(-EH(j,c(j)));
               end
               Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j,c(j))));
            elseif overdisp == 1
               if VH(j,c(j)) <= EH(j,c(j))
                  acomp = 0;
               else
                  acomp = 1;
                  a = EH(j,c(j))^2/(VH(j,c(j))-EH(j,c(j)));
               end
               if step == 1
                  Lambda(i,k) = Lambda(i,k)*EH(j,c(j));
               elseif step == -1
                  if acomp == 0
                     Lambda(i,k) = Lambda(i,k)*EH(j,c(j))*exp(-EH(j,c(j)));
                  elseif acomp == 1
                     Lambda(i,k) = Lambda(i,k)*EH(j,c(j))*(1+EH(j,c(j))/a)^(-a-1);
                  end
               end
               if acomp == 0
                  Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j,c(j))));
               elseif acomp == 1
                  Lambda(i,k) = Lambda(i,k)/(1-(1+EH(j,c(j))/a)^(-a));
               end
            end
            sratio = s(g(k)+1)/s(g(i)+1);
            if sratio==1
               Lambda(i,k) = Lambda(i,k)/N;
            else 
               Lambda(i,k) = Lambda(i,k)*(1-1/sratio)/(1-1/sratio^N);
            end
            if g(k)<g(i)
               Lambda(i,k) = Lambda(i,k)*gamma;
            end
%            i
%            k
%            Lambda(i,k)
         end
      end
      if tunnel == 1
         for step = [-2 2]
            d = c;
            d(j) = c(j)+step;
            if (max(d)<=C) & (min(d)>=1)
               k = lexorder(d,C,m);
               Lambda(i,k) = 2*N*mu*(m*L-1)^(-1)*M2(c(j),d(j),j)/M(j,c(j));
               if overdisp == 0
                  if step == 2
                     Lambda(i,k) = Lambda(i,k)*EH(j,c(j)); 
                  elseif step == -2
                     Lambda(i,k) = Lambda(i,k)*EH(j,c(j))*exp(-(EH(j,c(j))+EH(j,c(j)-1)));
                  end
                  Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j,c(j))));
               elseif overdisp == 1
                  if VH(j,c(j)) <= EH(j,c(j))
                     acomp = 0;
                  else
                     acomp = 1;
                     a = EH(j,c(j))^2/(VH(j,c(j))-EH(j,c(j)));
                  end
                  if step == 2
                     Lambda(i,k) = Lambda(i,k)*EH(j,c(j));
                  elseif step == -2
                     if acomp == 0
                        Lambda(i,k) = Lambda(i,k)*EH(j,c(j))*exp(-EH(j,c(j)));
                     elseif acomp == 1
                        Lambda(i,k) = Lambda(i,k)*EH(j,c(j))*(1+EH(j,c(j))/a)^(-a-1);
                     end
                     if EH(j,c(j)-1) <= VH(j,c(j)-1)
                        Lambda(i,k) = Lambda(i,k)*exp(-EH(j,c(j)-1));
                     else
                        a2 =  EH(j,c(j)-1)^2/(VH(j,c(j)-1)-EH(j,c(j)-1));
                        Lambda(i,k) = Lambda(i,k)*(1+EH(j,c(j)-1)/a2)^(-a2);
                     end
                  end
                  if acomp == 0
                     Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j,c(j))));
                  elseif acomp == 1
                     Lambda(i,k) = Lambda(i,k)/(1-(1+EH(j,c(j))/a)^(-a));
                  end
               end
               if s(g(k)+1)==s(g(i)+1)
                  beta = 1/N;
               else 
                  sratio = s(g(k)+1)/s(g(i)+1);
                  beta = (1-1/sratio)/(1-1/sratio^N);
               end
               Lambda(i,k) = Lambda(i,k)*beta;
               if g(k)<g(i)
                  Lambda(i,k) = Lambda(i,k)*gamma;
               end
               dpr = (d+c)/2;
               kpr = lexorder(dpr,C,m);
               sratiopr = s(g(kpr)+1)/s(g(i)+1);
               if sratiopr==1
                  betapr = 1/N;
               else
                  betapr = (1-1/sratiopr)/(1-1/sratiopr^N);
               end
               psi = (1-1/(3*(m*L-1)))*betapr + beta/(3*(m*L-1));
               r = (sratiopr-1 + sqrt((sratiopr-1)^2+2*(1+sratiopr)*sratiopr*(m*L-1)*mu*psi))/((1+sratiopr)*psi);
               Lambda(i,k) = Lambda(i,k)*r;
%            i
%            k
%            Lambda(i,k)
            end
         end
      end
   end
   if tunnel == 1
      for j1=1:(m-1)
         for j2=(j1+1):m
            for step1 = [-1 1]
               for step2 = [-1 1]
                  d = c;
                  d(j1) = c(j1)+step1;
                  d(j2) = c(j2)+step2;
                  if (max(d)<=C) & (min(d)>=1)
                     k = lexorder(d,C,m);
                     Lambda(i,k) = N*mu*(m*L-1)^(-1);
                     Lambda(i,k) = Lambda(i,k)*M1(c(j1),d(j1),j1)/M(j1,c(j1));
                     Lambda(i,k) = Lambda(i,k)*M1(c(j2),d(j2),j2)/M(j2,c(j2));
                     if overdisp == 0
                        if step1 == 1
                           Lambda(i,k) = Lambda(i,k)*EH(j1,c(j1));
                        elseif step1 == -1
                           Lambda(i,k) = Lambda(i,k)*EH(j1,c(j1))*exp(-EH(j1,c(j1)));
                        end
                        Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j1,c(j1))));
                        if step2 == 1
                           Lambda(i,k) = Lambda(i,k)*EH(j2,c(j2));
                        elseif step2 == -1
                           Lambda(i,k) = Lambda(i,k)*EH(j2,c(j2))*exp(-EH(j2,c(j2)));
                        end
                        Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j2,c(j2))));
                     elseif overdisp == 1
                        if VH(j1,c(j1)) <= EH(j1,c(j1))
                           acomp = 0;
                        else
                           acomp = 1;
                           a = EH(j1,c(j1))^2/(VH(j1,c(j1))-EH(j1,c(j1)));
                        end
                        if step == 1
                           Lambda(i,k) = Lambda(i,k)*EH(j1,c(j1));
                        elseif step == -1
                           if acomp == 0
                              Lambda(i,k) = Lambda(i,k)*EH(j1,c(j1))*exp(-EH(j1,c(j1)));
                           elseif acomp == 1
                              Lambda(i,k) = Lambda(i,k)*EH(j1,c(j1))*(1+EH(j1,c(j1))/a)^(-a-1);
                           end
                        end
                        if acomp == 0
                           Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j1,c(j1))));
                        elseif acomp == 1
                           Lambda(i,k) = Lambda(i,k)/(1-(1+EH(j1,c(j1))/a)^(-a));
                        end
                        if VH(j2,c(j2)) <= EH(j2,c(j2))
                           acomp = 0;
                        else
                           acomp = 1;
                           a = EH(j2,c(j2))^2/(VH(j2,c(j2))-EH(j2,c(j2)));
                        end
                        if step == 1
                           Lambda(i,k) = Lambda(i,k)*EH(j2,c(j2));
                        elseif step == -1
                           if acomp == 0
                              Lambda(i,k) = Lambda(i,k)*EH(j2,c(j2))*exp(-EH(j2,c(j2)));
                           elseif acomp == 1
                              Lambda(i,k) = Lambda(i,k)*EH(j2,c(j2))*(1+EH(j2,c(j2))/a)^(-a-1);
                           end
                        end
                        if acomp == 0
                           Lambda(i,k) = Lambda(i,k)/(1-exp(-EH(j2,c(j2))));
                        elseif acomp == 1
                           Lambda(i,k) = Lambda(i,k)/(1-(1+EH(j2,c(j2))/a)^(-a));
                        end
                     end
                     if s(g(k)+1)==s(g(i)+1)
                        beta = 1/N;
                     else 
                        sratio = s(g(k)+1)/s(g(i)+1);
                        beta = (1-1/sratio)/(1-1/sratio^N);
                     end
                     Lambda(i,k) = Lambda(i,k)*beta;
                     if g(k)<g(i)
                        Lambda(i,k) = Lambda(i,k)*gamma;
                     end
                     d1=c;
                     d1(j1)=d(j1);
                     k1 = lexorder(d1,C,m);
                     sratio1 = s(g(k1)+1)/s(g(i)+1);
                     if sratio1 == 1
                        beta1 = 1/N;
                     else 
                        beta1 = (1-1/sratio1)/(1-1/sratio1^N);
                     end
                     psi1 = (1-1/(3*(m*L-1)))*beta1 + beta/(3*(m*L-1));
                     r1 = (sratio1-1 + sqrt((sratio1-1)^2+2*(1+sratio1)*sratio1*(m*L-1)*mu*psi1))/((1+sratio1)*psi1);
                     d2=c;
                     d2(j2)=d(j2);
                     k2 = lexorder(d2,C,m);
                     sratio2 = s(g(k2)+1)/s(g(i)+1);
                     if sratio2 == 1
                        beta2 = 1/N;
                     else 
                        beta2 = (1-1/sratio2)/(1-1/sratio2^N);
                     end
                     psi2 = (1-1/(3*(m*L-1)))*beta2 + beta/(3*(m*L-1));
                     r2 = (sratio2-1 + sqrt((sratio2-1)^2+2*(1+sratio2)*sratio2*(m*L-1)*mu*psi2))/((1+sratio2)*psi2);
                     Lambda(i,k) = Lambda(i,k)*(r1+r2);
 %                 i
 %                 k
 %                 Lambda(i,k)
                  end
               end
            end
         end
      end
   end
   Lambda(i,i) = - sum(Lambda(i,:));
end
%Lambda
Lambda0 = Lambda(Inonabs,Inonabs);
lambda1 = Lambda(Inonabs,Iabs)*ones(length(Iabs),1);

% The wating time is defined as the time when the array component proccess C(t) reaches an 
% absorbing state, i.e. a state c=[c(1),...,c(m)] such that g(c)=m. Compute expected value, standard 
% deviation, density function, distribution function and quantiles of the waiting time distribution, 
% using the fact that C(t) is assumed to be a continuous time Markov process and hence T has a 
% phase-type distribution.
onevec = ones(length(Inonabs),1);
ET = - kappa0'*inv(Lambda0)*onevec;
DT = sqrt(2*kappa0'*inv(Lambda0)^2*onevec-ET^2);
if length(tvec)==0
      fT = [];
      FT = [];
else
   for i=1:length(tvec)
      fT(i,1) = kappa0'*expm(Lambda0*tvec(i))*lambda1;
      FT(i,1) = 1 - kappa0'*expm(Lambda0*tvec(i))*onevec;
   end
end
if length(qvec)==0
   qT = [];
else
   if length(Inonabs)==1 % Waiting time has an exponential distribution if there is only one non-absorbing state.
      qT = (log(min([ones(1,length(qvec));(1-qvec)/kappa0'])))'/Lambda0;
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
               F = 1 - kappa0'*expm(Lambda0*qTmid)*onevec;
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

% Adjust output parameter
g = g';
