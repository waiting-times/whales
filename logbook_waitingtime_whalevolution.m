% Numerical illustrations for paper 'Temporal conttraints for the evolution of the whale locomotory
% system" by Bechly, Gauger, Sternberg and Hössjer.
% 
% The calculations are based on the Matlab programs waitingtime_regseq_2.m and waitingtime_regseq_iter_2.m 
% used in the paper Hössjer, Bechly and Gauger (2020), which will be submitted to Journal of Theoretical 
% Biology. In order to include generation time Deltat as input parameter, and modify the program
% to handle large K and m (see below) we also implemented the Matlab probrams 
% waitingtime_regseq_whalevolution and waitingtime_regseq_whalevolution_iter. 

% Default parameter values (order of appearance as for the input parameters of 
% waitingtime_regseq_whalevolution.m and waitingtime_regseq_whalevolution_iter.m):
% Generation time Deltat = 5 (years) 
% Haploid population size N = 2*50,000=100,000
% Length of regulatory sequence of each gene L = 1000 or 5000
% Number of possible genes M = m
% Number of genes m varies between 1 and 500
% Mutation rate mu = 2.2*10^(-9) (default value) or 10^(-8)
% Back mutation probability gamma = 0 or 1
% Binding site length W = 8,...,15 (W=10 default)
% Fitness type fitntype = 'Final'
% Selection coefficients s = ones(1,m+1) (neutral model) or s = [1 1 ... 1 10] (target-selected model). 
% Number of binding sites per gene K = 3, 10, 20 or 30 (K=3 default value) 
% Maximal number of mismatches deltamax = 0,1,2,3,4 (deltamax = 2 default)
% Number of distance intervals per gene C = 2 for neutral model, C = 3 (or C = 2, for large m) for target selected model.
% TA = 'arbitrary' (not an input parameter, default in waitingtime_regseq_whalevolution.m)
% tvec = tmax = 1.2*10^6 years is the length of the time window)
% qvec = 0.5
% tunnel = 0 for neutral model, = 1 for target selected model 

% Preparation: Ration between exptected value and median of an exponential distribuiton
% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
>> -1/log(0.5)

ans =

    1.4427

% Initial run with waitingtime_regseq_2
% -------------------------------------
% -------------------------------------
>> Deltat = 5; 
>> W = 8;
>> b = [2*ones(1,W); 3*ones(1,W)]; 
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_2(10^5,1000,1,2.2*10^(-9),0,8,[1 1],3,[0 1],0,'arbitrary',2.4*10^5,[0.1 0.5 0.9],b,1)
>> ET = ET/Deltat
>> DT = DT/Deltat
>> qT = qT/Deltat
>> FT
% This run was interrupted, because it takes too long with W=8 and K=3.

% A number of runs with waitingtime_regseq_whalevolution_iter.m
% -------------------------------------------------------------
% -------------------------------------------------------------

% 1) W = 10, K = 3, deltamax = 2, m varies
% -----------------------------------------
% ----------------------------------------

% 1a) Neutral model, no back mutations
% ------------------------------------

% 1aI) C = 2, tunnel = 0 (numerical results of FT unstable for large m due to subtraction of almost
% equally large numbers)
% -------------------------------------------------------------------------------------------------
>> m = [1 2 3 4 5 10 20 50 100 200];
>> gamma = 0;
>> s = 1;
>> C = 2;
>> tunnel = 0;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+009 *

    0.0719
    0.1333
    0.1863
    0.2325
    0.2731
    0.4203
    0.5843
    0.8073
    0.9775
    1.1484


DT =

  1.0e+008 *

    1.7435
    2.2533
    2.5408
    2.7199
    2.8370
    3.0574
    3.1252
    3.1539
    3.1635
    3.1683


fT =

  1.0e-008 *

    0.1168    0.1661    0.1772    0.1679    0.1492    0.0542    0.0036    0.0000    0.0000    0.0000


FT =

    0.7109    0.5054    0.3593    0.2554    0.1816    0.0330    0.0011    0.0000    0.0000         0


qT =

  1.0e+009 *

         0         0    0.0847    0.1489    0.2000    0.3630    0.5303    0.7544    0.9250    1.0961

>> FT(7)

ans =

    0.0011

>> FT(8)

ans =

  3.8980e-008

>> FT(9)

ans =

  1.3323e-015

>> FT(10)

ans =

     0
     
>> qT(1)

ans =

     0

>> qT(2)

ans =

     0

>> qT(3)

ans =

  8.4673e+007
  
% Adjust probability FT(9) for m=100 by intgrating density function 
>> [ET,DT,fT,FT,qT,kappa,lambda,FT2] = waitingtime_regseq_whaleevolution(5,100000,1000,100,100,2.2*10^(-9),0,10,'Final',1,3,2,2,1.2*10^6,0.5,tunnel,50);
>> FT % This is FT(9)

FT =

  1.3323e-015

>> FT2 % This is an improved calculation of FT(9) based on numerical integration

FT2 =

  1.5194e-015
     
% 1aII) C = 3, tunnel = 1 (makes hardly any difference to assume C=3 and tunnel = 1 for neutral model)
% ----------------------------------------------------------------------------------------------------
>> m = [1 2 3 4 5 10 20];
>> gamma = 0;
>> s = 1;
>> C = 3;
>> tunnel = 1;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    0.7187
    1.3330
    1.8632
    2.3249
    2.7308
    4.2032
    5.8430


DT =

  1.0e+008 *

    1.7435
    2.2533
    2.5409
    2.7199
    2.8370
    3.0575
    3.1252


fT =

  1.0e-008 *

    0.1168    0.1661    0.1772    0.1679    0.1492    0.0542    0.0036


FT =

    0.7109    0.5054    0.3593    0.2554    0.1816    0.0330    0.0011


qT =

  1.0e+008 *

         0         0    0.8467    1.4894    1.9998    3.6304    5.3028

>> qT(1)

ans =

     0

>> qT(2)

ans =

     0

% 1b) Neutral model, back mutations
% ---------------------------------

% 1bI) C = 2 and small values of m
% --------------------------------
>> m = [1 2 3 4 5 10 20];
>> gamma=1;
>> s = 1;
>> C = 2;
>> tunnel = 0;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+010 *

    0.0072
    0.0158
    0.0264
    0.0393
    0.0553
    0.2259
    3.4007


DT =

  1.0e+010 *

    0.0174
    0.0280
    0.0392
    0.0522
    0.0681
    0.2358
    3.4053


fT =

  1.0e-008 *

    0.1168    0.1658    0.1765    0.1669    0.1481    0.0532    0.0034


FT =

    0.7109    0.5054    0.3593    0.2554    0.1816    0.0330    0.0011


qT =

  1.0e+010 *

         0         0    0.0096    0.0198    0.0319    0.1531    2.3557

>> ET(1)

ans =

  7.1869e+007

>> FT(7)

ans =

    0.0011

>> qT(1)

ans =

     0

>> qT(2)

ans =

     0

>> qT(3)

ans =

  9.5600e+007
  
% 1bII) C = 2 and large m (numerical results of ET and DT for m =100 and 200 unstable due to inversion
% of matrix Lambda0, numerical results for FT unstable due to cancellatin, i.e. subtraction of almost
% equally large numbers, numerical results of qT unstable for m=200, since qT does not increasase
% as much as it should compared to m=100)
% --------------------------------------------------------------------------------------------------------
>> m = [50 100 200];
>> gamma=1;
>> s=1;
>> C=2;
>> tunnel=0;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+023 *

    0.0000
    0.0510
    3.6388


DT =

  1.0e+023 *

    0.0000
    0.0510
    3.6388


fT =

  1.0e-014 *

    0.2908    0.0000    0.0000


FT =

  1.0e-007 *

    0.3880    0.0000         0


qT =

  1.0e+025 *

    0.0000    0.0003    1.1884

>> ET(1)

ans =

  3.6144e+014

>> ET(2)

ans =

  5.0955e+021

>> FT(2)

ans =

  1.4433e-015

>> FT(3)

ans =

     0

>> qT(1)

ans =

  2.5051e+014

>> qT(2)

ans =

  3.2976e+021

>> ET(1)/qT(1)

ans =

    1.4428

>> ET(2)/qT(2)

ans =

    1.5452

>> ET(3)/qT(3)

ans =

    0.0306

% Corrected values of ET(2) and ET(3), based on assumption that T has an exponential distribution.
>> ETcorr2 = qT(2)/log(2)

ETcorr2 =

  4.7575e+021

>> ETcorr3 = qT(3)/log(2)

ETcorr3 =

  1.7144e+025
  
% Compute probability P(T=0)
>> m=100;
>> gamma=1;
>> s=1;
>> C=2;
>> tunnel=0;
>> [ET,DT,fT,FT,qT,kappa] = waitingtime_regseq_whaleevolution(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)
>> kappa(101) % This is P(T=0)

ans =

  1.2466e-015

>> kappa(100) % This is probability of one missing binding site at time t=0.

ans =

  5.1040e-014

>> kappa(99) % This is probability of two missing binding sites at time t=0.

ans =

  1.0344e-012
  
% Adjust probability FT(2) for m=100 by integrating density function
>> [ET,DT,fT,FT,qT,kappa,lambda,FT2] = waitingtime_regseq_whaleevolution(5,100000,1000,100,100,2.2*10^(-9),1,10,'Final',1,3,2,2,1.2*10^6,0.5,tunnel,50);
>> FT  % This is FT(2)

FT =

  1.4433e-015

>> FT2 % This is FT(2) computed more accurately by numerical integration.

FT2 =

  1.4936e-015

>> m=100;
>> gamma=1;
>> s=1;
>> C=2;
>> tunnel=0;
>> tmax = 1.2*10^6;
>> [ET,DT,fT,FT,qT,kappa] = waitingtime_regseq_whaleevolution(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,(1:20)*tmax/20,0.5,tunnel);
>> FT2appr = kappa(101) + tmax*sum(fT)/20

FT2appr =

  1.4936e-015

% Verify that m=200 gives unreliable values of qT, by computing qT for smaller values of m and 
% noticing that qT is not an increasing function of m.
>> m=150;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel);
>> qT

qT =

  3.7101e+022

>> m=160;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel);
>> qT 

qT =

  1.2788e+025
  
>> m=170;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel);
>> qT

qT =

  2.2146e+025 
  
>> m=200;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel);
>> qT

qT =

  1.1884e+025

% 1bIII) C = 3 and small m (values essentially identical to those for C=2 in 1bI) 
% -------------------------------------------------------------------------------
>> m = [1 2 3 4 5 10 20];
>> gamma=1;
>> s=1;
>> C=3;
>> tunnel=1;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+010 *

    0.0072
    0.0158
    0.0264
    0.0393
    0.0553
    0.2259
    3.4007


DT =

  1.0e+010 *

    0.0174
    0.0280
    0.0392
    0.0522
    0.0681
    0.2358
    3.4053


fT =

  1.0e-008 *

    0.1168    0.1658    0.1765    0.1669    0.1481    0.0532    0.0034


FT =

    0.7109    0.5054    0.3593    0.2554    0.1816    0.0330    0.0011


qT =

  1.0e+010 *

         0         0    0.0096    0.0198    0.0319    0.1531    2.3555

ET(1)

ans =

  7.1870e+007

FT(7)

ans =

    0.0011

qT(1)

ans =

     0

qT(2)

ans =

     0

qT(3)

ans =

  9.5610e+007
  
% 1c) Target selected model, no back mutations
% --------------------------------------------

% 1cI) C = 3, tunnel = 0
% ----------------------
>> m = [1 2 3 4 5 10 20];
>> gamma=0;
>> s = 10;
>> C = 3;
>> tunnel = 0;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    0.0000
    0.1044
    0.2728
    0.4779
    0.7016
    1.8091
    3.3715


DT =

  1.0e+008 *

    0.0027
    0.4974
    0.7961
    1.0349
    1.2253
    1.7159
    1.9049


fT =

  1.0e-008 *

    0.0001    0.0677    0.1442    0.2050    0.2428    0.1984    0.0276


FT =

    1.0000    0.9164    0.7976    0.6708    0.5507    0.1670    0.0099


qT =

  1.0e+008 *

         0         0    0.0000    0.0000    0.0001    1.4409    3.1141

>> ET(1)

ans =

  1.8502e+003
  
>> qT(1)

ans =

     0

>> qT(2)

ans =

     0

>> qT(3)

ans =

  1.0828e+003

>> qT(4)

ans =

  2.4791e+003

>> qT(5)

ans =

  5.6080e+003
  
% Additional run with qvec = [0.4 0.5 0.6]
>> qT =

  1.0e+008 *

         0         0    0.0000    0.0000    0.0000    1.0227    2.6742
         0         0    0.0000    0.0000    0.0001    1.4409    3.1141
         0    0.0000    0.0000    0.0000    0.2221    1.8935    3.5859

>> qT(1,5)

ans =

  2.5070e+003

>> qT(2,5)

ans =

  5.6080e+003

>> qT(3,5)

ans =

  2.2206e+007

% 1cII) C = 2 (for comparison with C = 3), tunnel = 0 (values only marginally larger than in 1cI)
% -----------------------------------------------------------------------------------------------
>> m = [1 2 3 4 5 10 20];
>> gamma=0;
>> s = 10;
>> C = 2;
>> tunnel = 0;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    0.0000
    0.1044
    0.2727
    0.4778
    0.7015
    1.8091
    3.3714


DT =

  1.0e+008 *

    0.0000
    0.4974
    0.7961
    1.0348
    1.2253
    1.7158
    1.9048


fT =

  1.0e-008 *

    0.0000    0.0676    0.1441    0.2049    0.2427    0.1983    0.0276


FT =

    1.0000    0.9164    0.7976    0.6709    0.5508    0.1671    0.0099


qT =

  1.0e+008 *

         0         0    0.0000    0.0000    0.0001    1.4408    3.1140  


% 1cIII) C = 3, tunnel = 1 (the values are an order of magnitute larger with stochastic tunneling)
% -----------------------------------------------------------------------------------------------
>> m = [1 2 3 4 5 10 20];
>> gamma = 0;
>> s = 10;
>> C = 3;
>> tunnel = 1;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    0.0000
    0.1044
    0.2728
    0.4779
    0.7016
    1.8091
    3.3715


DT =

  1.0e+008 *

    0.0027
    0.4974
    0.7961
    1.0349
    1.2253
    1.7159
    1.9049


fT =

  1.0e-008 *

    0.0001    0.0677    0.1442    0.2050    0.2428    0.1984    0.0276


FT =

    1.0000    0.9164    0.7976    0.6708    0.5507    0.1670    0.0099


qT =

  1.0e+008 *

         0         0    0.0000    0.0000    0.0001    1.4409    3.1141

>> ET(1)

ans =

  1.8502e+003

>> FT(7)

ans =

    0.0099

>> qT(1)

ans =

     0

>> qT(2)

ans =

     0

>> qT(3)

ans =

  1.0828e+003

>> qT(4)

ans =

  2.4791e+003

>> qT(5)

ans =

  5.6080e+003

>> ET(7)/qT(7)

ans =

    1.4528

ET(6)/qT(6)

ans =

    1.7454

ET(5)/qT(5)

ans =

  1.8014e+004
         
% 1d) Target selected model, back mutations
% -----------------------------------------
>> m = [1 2 3 4 5 10 20];
>> gamma = 1;
>> s = 10;
>> C = 3;
>> tunnel = 1;
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+009 *

    0.0000
    0.0104
    0.0307
    0.0608
    0.1010
    0.5087
    5.1572


DT =

  1.0e+009 *

    0.0003
    0.0497
    0.0927
    0.1405
    0.1943
    0.6261
    5.2332


fT =

  1.0e-008 *

    0.0001    0.0677    0.1439    0.2042    0.2414    0.1953    0.0267


FT =

    1.0000    0.9164    0.7976    0.6708    0.5507    0.1670    0.0099


qT =

  1.0e+009 *

         0         0    0.0000    0.0000    0.0000    0.2915    3.5499

>> ET(1)

ans =

  1.8502e+003

>> FT(7)

ans =

    0.0099

>> qT(1)

ans =

     0

>> qT(2)

ans =

     0

>> qT(3)

ans =

  1.0828e+003

>> qT(4)

ans =

  2.4791e+003

>> qT(5)

ans =

  5.6083e+003

% 2) K = 3, deltmax = 2, m = 5, W varies
% ---------------------------------------
% ---------------------------------------

% 2a) Neutral model, no back mutations
% ------------------------------------
>> m = 5;
>> gamma = 0;
>> s = 1;
>> C = 3;
>> tunnel = 1;
>> W = [8 10 12 15];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,W,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+011 *

    0.0000
    0.0027
    0.0489
    1.6551


DT =

  1.0e+010 *

    0.0000
    0.0284
    0.2727
    8.7805


fT =

  1.0e-008 *

    0.0000    0.1492    0.0000    0.0000


FT =

    1.0000    0.1816    0.0000    0.0000


qT =

  1.0e+011 *

         0    0.0020    0.0435    1.4817

>> ET(1)

ans =

  562.3987

>> ET(2)

ans =

  2.7308e+008

>> FT(3)

ans =

  1.3436e-005

>> FT(4)

ans =

  1.5554e-013

>> qT(1)

ans =

     0

>> qT(2)

ans =

  1.9998e+008

>> ET(4)/qT(4)

ans =

    1.1170

>> ET(3)/qT(3)

ans =

    1.1238

>> ET(2)/qT(2)

ans =

    1.3656


% 2b) Neutral model, back mutations
% ---------------------------------
>> m=5;
>> gamma=1;
>> s=1;
>> C=3;
>> tunnel=1;
>> W = [8 10 12 15];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,W,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+020 *

    0.0000
    0.0000
    0.0000
    2.4793


DT =

  1.0e+020 *

    0.0000
    0.0000
    0.0000
    2.4793


fT =

  1.0e-008 *

    0.0000    0.1481    0.0000    0.0000


FT =

    1.0000    0.1816    0.0000    0.0000


qT =

  1.0e+020 *

         0    0.0000    0.0000    1.7181

>> ET(1)

ans =

  562.4025

>> ET(2)

ans =

  5.5321e+008

>> ET(3)

ans =

  4.0671e+012

>> FT(3)

ans =

  1.3433e-005

>> FT(4)

ans =

  1.5576e-013

>> qT(1)

ans =

     0

>> qT(2)

ans =

  3.1895e+008

>> qT(3)

ans =

  2.8191e+012

>> ET(4)/qT(4)

ans =

    1.4431

>> ET(3)/qT(3)

ans =

    1.4427

>> ET(2)/qT(2)

ans =

    1.7345

% Refined check if FT(4) is numerically unreliable du to cancellation when terms are subtracted.
>> [ET,DT,fT,FT,qT,kappa,lambda,FT2] = waitingtime_regseq_whaleevolution(5,100000,1000,5,5,2.2*10^(-9),1,15,'Final',1,3,2,3,1.2*10^6,0.5,tunnel,50);
>> FT

FT =

  1.5576e-013

>>FT2

FT2 =

  1.5555e-013

% 2c) Target selected model, no back mutations
% --------------------------------------------
>> m=5;
>> gamma=0;
>> s=10;
>> C=3;
>> tunnel=1;
>> W = [8 10 12 15];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,W,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+010 *

    0.0000
    0.0070
    0.2725
    9.8355


DT =

  1.0e+010 *

    0.0000
    0.0123
    0.1546
    4.9733


fT =

  1.0e-008 *

    0.0000    0.2428    0.0007    0.0000


FT =

    1.0000    0.5507    0.0004    0.0000


qT =

  1.0e+010 *

         0    0.0000    0.2448    8.9439

>> ET(1)

ans =

    0.0082

>> ET(2)

ans =

  7.0159e+007

>> FT(3)

ans =

  3.8228e-004

>> FT(4)

ans =

  9.6098e-012

>> qT(1)

ans =

     0

>> qT(2)

ans =

  5.6080e+003

>> ET(4)/qT(4)

ans =

    1.0997

>> ET(3)/qT(3)

ans =

    1.1132

>> ET(2)/qT(2)

ans =

  1.2510e+004
  
% 2d) Target selected model, back mutations
% -----------------------------------------
>> m=5;
>> gamma=1;
>> s=10;
>> C=3;
>> tunnel=1;
>> W = [8 10 12 15];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,W,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+018 *

    0.0000
    0.0000
    0.0000
    4.1657


DT =

  1.0e+018 *

    0.0000
    0.0000
    0.0000
    4.1657


fT =

  1.0e-008 *

    0.0000    0.2414    0.0006    0.0000


FT =

    1.0000    0.5507    0.0004    0.0000


qT =

  1.0e+018 *

         0    0.0000    0.0000    2.8782

>> ET(1)

ans =

    0.0082

>> ET(2)

ans =

  1.0103e+008

>> ET(3)

ans =

  1.7526e+011

>> ET(4)

ans =

  4.1657e+018

>> FT(3)

ans =

  3.8213e-004

>> FT(4)

ans =

  9.6021e-012

>> qT(1)

ans =

     0

>> qT(2)

ans =

  5.6083e+003

>> qT(3)

ans =

  1.2145e+011

>> ET(2)/qT(2)

ans =

  1.8014e+004

>> ET(3)/qT(3)

ans =

    1.4430

>> ET(4)/qT(4)

ans =

    1.4473


% 3) W = 10, deltamax = 2, m = 5, K varies
% ----------------------------------------
% ----------------------------------------

% 3a) Neutral model, no back mutations
% ------------------------------------
>> m=5;
>> gamma=0;
>> s=1;
>> C=3;
>> tunnel=1;
>> K = [3 10 20 30];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,K,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    2.7308
    0.0593
    0.0005
    0.0000


DT =

  1.0e+008 *

    2.8370
    0.2919
    0.0190
    0.0016


fT =

  1.0e-008 *

    0.1492    0.1009    0.0034    0.0001


FT =

    0.1816    0.9226    0.9987    1.0000


qT =

  1.0e+008 *

    1.9998         0         0         0

>> ET(2)

ans =

  5.9281e+006

>> ET(3)

ans =

  4.8893e+004

>>  ET(4)

ans =

  529.3181
  
>> qT(2)

ans =

     0

>> qT(3)

ans =

     0

>> qT(4)

ans =

     0
     
>> ET(1)/qT(1)

ans =

    1.3656

% 3b) Neutral model, back mutations
% ---------------------------------
>> m=5;
>> gamma=1;
>> s=1;
>> C=3;
>> tunnel=1;
>> K = [3 10 20 30];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,K,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    5.5321
    0.0613
    0.0005
    0.0000


DT =

  1.0e+008 *

    6.8066
    0.3041
    0.0191
    0.0016


fT =

  1.0e-008 *

    0.1481    0.1008    0.0034    0.0001


FT =

    0.1816    0.9226    0.9987    1.0000


qT =

  1.0e+008 *

    3.1895         0         0         0

>> ET(2)

ans =

  6.1254e+006

>> ET(3)

ans =

  4.8919e+004

>> ET(4)

ans =

  529.3226

>> qT(2)

ans =

     0

>> qT(3)

ans =

     0

>> qT(4)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.7345

% 3c) Target selected model, no back mutations
% --------------------------------------------
>> m=5;
>> gamma=0;
>> s=10;
>> C=3;
>> tunnel=1;
>> K = [3 10 20 30]; 
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,K,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+007 *

    7.0159
    0.0096
    0.0000
    0.0000


DT =

  1.0e+008 *

    1.2253
    0.0267
    0.0002
    0.0000


fT =

  1.0e-008 *

    0.2428    0.0066    0.0000    0.0000


FT =

    0.5507    0.9975    1.0000    1.0000


qT =

  1.0e+003 *

    5.6080         0         0         0

>> ET(2)

ans =

  9.5779e+004

>> ET(3)

ans =

   13.4287

>> ET(4)

ans =

    0.0081

>> qT(2)

ans =

     0

>> qT(3)

ans =

     0

>> qT(4)

ans =

     0

>> ET(1)/qT(1)

ans =

  1.2510e+004

% 3d) Target selected model, back mutations
% -----------------------------------------
>> m=5;
>> gamma=1;
>> s=10;
>> C=3;
>> tunnel=1;
>> K = [3 10 20 30];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,K,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    1.0103
    0.0010
    0.0000
    0.0000


DT =

  1.0e+008 *

    1.9432
    0.0273
    0.0002
    0.0000


fT =

  1.0e-008 *

    0.2414    0.0065    0.0000    0.0000


FT =

    0.5507    0.9975    1.0000    1.0000


qT =

  1.0e+003 *

    5.6083         0         0         0

>> ET(2)

ans =

  9.7365e+004

>> ET(3)

ans =

   13.4324

>> ET(4)

ans =

    0.0081

>> qT(2)

ans =

     0

>> qT(3)

ans =

     0

>> qT(4)

ans =

     0

>> ET(1)/qT(1)

ans =

  1.8014e+004


% 4) W = 10, K = 3, m = 5, deltamax varies
% ----------------------------------------
% ----------------------------------------

% 4a) Neutral model, no back mutations
% ------------------------------------
>> m=5;
>> gamma=0;
>> s=1;
>> C=3;
>> tunnel=1;
>> deltamax = [0 1 2 3 4];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,deltamax,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+011 *

    1.8903
    0.0670
    0.0027
    0.0000
    0.0000


DT =

  1.0e+011 *

    1.0028
    0.0369
    0.0028
    0.0000
    0.0000


fT =

  1.0e-008 *

    0.0000    0.0000    0.1492    0.0004    0.0000


FT =

    0.0000    0.0000    0.1816    0.9999    1.0000


qT =

  1.0e+011 *

    1.6923    0.0597    0.0020         0         0

>> ET(3)

ans =

  2.7308e+008

>> ET(4)

ans =

  5.2583e+003

>> ET(5)

ans =

  1.3258e-018

>> FT(1)

ans =

  1.8663e-013

>> FT(2)

ans =

  4.3110e-006

>> qT(3)

ans =

  1.9998e+008

>> qT(4)

ans =

     0

>> qT(5)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.1171

>> ET(2)/qT(2)

ans =

    1.1222

>> ET(3)/qT(3)

ans =

    1.3656
    
% Check if FT(1) is numerically reliable due to cancellation. The numercial integration computation
% of FT(1) based on Simpson's formula (=FT2) agrees up to three digits. 
>> [ET,DT,fT,FT,qT,kappa,lambda,FT2] = waitingtime_regseq_whaleevolution(5,100000,1000,5,5,2.2*10^(-9),0,10,'Final',1,3,0,3,1.2*10^6,0.5,tunnel,50);
>> FT

FT =

  1.8663e-013

>> FT2

FT2 =

  1.8677e-013
  
% 4b) Neutral model, back mutations
% ---------------------------------
>> m=5;
>> gamma=1;
>> s=1;
>> C=3;
>> tunnel=1;
>> deltamax = [0 1 2 3 4];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,deltamax,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+020 *

    2.5204
    0.0000
    0.0000
    0.0000
    0.0000


DT =

  1.0e+020 *

    2.5204
    0.0000
    0.0000
    0.0000
    0.0000


fT =

  1.0e-008 *

    0.0000    0.0000    0.1481    0.0004    0.0000


FT =

    0.0000    0.0000    0.1816    0.9999    1.0000


qT =

  1.0e+020 *

    1.7470    0.0000    0.0000         0         0

>> ET(2)

ans =

  1.3360e+013

>> ET(3)

ans =

  5.5321e+008

>> ET(4)

ans =

  5.2586e+003

>> ET(5)

ans =

  1.3258e-018

>> FT(1)

ans =

  1.8652e-013

>> FT(2)

ans =

  4.3102e-006

>> qT(2)

ans =

  9.2607e+012

>> qT(3)

ans =

  3.1895e+008

>> qT(4)

ans =

     0

>> qT(5)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.4427

>> ET(2)/qT(2)

ans =

    1.4427

>> ET(3)/qT(3)

ans =

    1.7345

% 4c) Target selected model, no back mutations
% --------------------------------------------
>> m=5;
>> gamma=0;
>> s=10;
>> C=3;
>> tunnel=1;
>> deltamax = [0 1 2 3 4];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,deltamax,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+011 *

    1.0887
    0.0373
    0.0007
    0.0000
    0.0000


DT =

  1.0e+010 *

    5.6519
    0.2086
    0.0123
    0.0000
    0.0000


fT =

  1.0e-008 *

    0.0000    0.0003    0.2428    0.0000    0.0000


FT =

    0.0000    0.0002    0.5507    1.0000    1.0000


qT =

  1.0e+010 *

    9.8595    0.3350    0.0000         0         0

>> ET(3)

ans =

  7.0159e+007

>> ET(4)

ans =

    0.2149

>> ET(5)

ans =

  1.4732e-023

>> FT(1)

ans =

  2.6932e-011

>> FT(2)

ans =

  1.6401e-004

>> qT(3)

ans =

  5.6080e+003

>> qT(4)

ans =

     0

>> qT(5)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.1042

>> ET(2)/qT(2)

ans =

    1.1125

>> ET(3)/qT(3)

ans =

  1.2510e+004

% 4d) Target selected model, back mutations
% -----------------------------------------
>> m=5;
>> gamma=1;
>> s=10;
>> C=3;
>> tunnel=1;
>> deltamax = [0 1 2 3 4];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,1000,m,m,2.2*10^(-9),gamma,10,'Final',s,3,deltamax,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+018 *

    1.8227
    0.0000
    0.0000
    0.0000
    0.0000


DT =

  1.0e+018 *

    1.8227
    0.0000
    0.0000
    0.0000
    0.0000


fT =

  1.0e-008 *

    0.0000    0.0003    0.2414    0.0000    0.0000


FT =

    0.0000    0.0002    0.5507    1.0000    1.0000


qT =

  1.0e+018 *

    1.2641    0.0000    0.0000         0         0

>> ET(2)

ans =

  4.2624e+011

>> ET(3)

ans =

  1.0103e+008

>> ET(4)

ans =

    0.2149

>> ET(5)

ans =

  1.4732e-023

>> FT(1)

ans =

  2.6893e-011

>> FT(2)

ans =

  1.6394e-004

>> qT(2)

ans =

  2.9543e+011

>> qT(3)

ans =

  5.6083e+003

>> qT(4)

ans =

     0

>> qT(5)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.4418

>> ET(2)/qT(2)

ans =

    1.4428

>> ET(3)/qT(3)

ans =

  1.8014e+004

% 5) W = 10, K = 3, m = 5, deltamax=2, L varies
% ---------------------------------------------
% ---------------------------------------------

% 5a) Neutral model, no back mutations
% ------------------------------------
>> m=5;
>> gamma=0;
>> s=1;
>> C=3;
>> tunnel=1;
>> L=[1000 5000];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,L,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    2.7308
    0.0048


DT =

  1.0e+008 *

    2.8370
    0.0689


fT =

  1.0e-008 *

    0.1492    0.0195


FT =

    0.1816    0.9904


qT =

  1.0e+008 *

    1.9998         0

>> ET(2)

ans =

  4.8481e+005

>> qT(2)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.3656
    
% 5b) Neutral model, back mutations
% ---------------------------------
>> m=5;
>> gamma=1;
>> s=1;
>> C=3;
>> tunnel=1;
>> L=[1000 5000];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,L,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    5.5321
    0.0049


DT =

  1.0e+008 *

    6.8066
    0.0692


fT =

  1.0e-008 *

    0.1481    0.0195


FT =

    0.1816    0.9904


qT =

  1.0e+008 *

    3.1895         0

>> ET(2)

ans =

  4.8673e+005

>> qT(2)

ans =

     0

>> ET(1)/qT(1)

ans =

    1.7345
    
% 5c) Target-selected model, no back mutations
% --------------------------------------------
>> m=5;
>> gamma=0;
>> s=10;
>> C=3;
>> tunnel=1;
>> L=[1000 5000];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,L,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+007 *

    7.0159
    0.0001


DT =

  1.0e+008 *

    1.2253
    0.0022


fT =

  1.0e-008 *

    0.2428    0.0002


FT =

    0.5507    1.0000


qT =

  1.0e+003 *

    5.6080         0

>> ET(2)

ans =

  963.5710

>> qT(2)

ans =

     0;

>> ET(1)/qT(1)

ans =

  1.2510e+004

% 5d) Target-selected model, back mutations
% -----------------------------------------
>> m=5;
>> gamma=1;
>> s=10;
>> C=3;
>> tunnel=1;
>> L=[1000 5000];
>> [ET,DT,fT,FT,qT] = waitingtime_regseq_whaleevolution_iter(5,100000,L,m,m,2.2*10^(-9),gamma,10,'Final',s,3,2,C,1.2*10^6,0.5,tunnel)

ET =

  1.0e+008 *

    1.0103
    0.0000


DT =

  1.0e+008 *

    1.9432
    0.0022


fT =

  1.0e-008 *

    0.2414    0.0002


FT =

    0.5507    1.0000


qT =

  1.0e+003 *

    5.6083         0

>> ET(2)

ans =

  965.4918

>> qT(2)

ans =

     0

>> ET(1)/qT(1)

ans =

  1.8014e+004
  
  