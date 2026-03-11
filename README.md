# Temporal Constraints for the Evolution of the Whale Locomotory System

This repository contains MATLAB scripts used for numerical illustrations in the paper:

Bechly, G., Gauger, A. K., Hössjer, O., Nelson, P. A., von Sternberg, R., and Luskin, C.  
*Temporal Constraints on the Evolution of the Whale Locomotory System (Mammalia: Cetacea): A Waiting Times Analysis*.


The calculations are adapted from MATLAB programs originally developed for:

Hössjer, Bechly, and Gauger (2020), submitted to the Journal of Theoretical Biology.

The original MATLAB programs:

- waitingtime_regseq_2.m
- waitingtime_regseq_iter_2.m

were modified to:

- include generation time (Δt) as an input parameter
- handle larger values of K (binding sites per gene) and m (number of genes)

The resulting MATLAB programs were:

- waitingtime_regseq_whalevolution.m
- waitingtime_regseq_whalevolution_iter.m

The MATLAB code in this repository reproduces the numerical calculations used in the whale locomotory system study.

---

# Model Parameters

The simulations use the following parameters (in the same order as the inputs of the original MATLAB functions).

## Population and Time

Δt  
Generation time.  
Default: 5 years

N  
Haploid population size.  
Default: 100,000

tmax  
Length of the time window.  
Default: 1.2 × 10^6 years

---

## Gene and Regulatory Sequence

L  
Length of regulatory sequence per gene.  
Values used: 1000 or 5000

M  
Number of possible genes.  
Defined as m.

m  
Number of genes.  
Values between 1 and 500.

---

## Mutation Parameters

μ  
Mutation rate.  
Default: 2.2 × 10^-9  
Alternative value: 1 × 10^-8

γ  
Back mutation probability.  
Values: 0 or 1

---

## Binding Site Model

W  
Binding site length.  
Values between 8 and 15.  
Default: 10

K  
Number of binding sites per gene.  
Values used: 3, 10, 20, or 30.  
Default: 3

Δmax (deltamax)  
Maximum number of mismatches.  
Values between 0 and 4.  
Default: 2

---

## Fitness Models

Neutral Model

fitntype = "Final"

Selection coefficients:
s = [1, 1, ..., 1]

Number of distance intervals per gene:
C = 2

tunnel = 0

---

Target-Selected Model

Selection coefficients:
s = [1, 1, ..., 1, 10]

Number of distance intervals per gene:
C = 3

For large m, the model may use:
C = 2

tunnel = 1

---

## Other Parameters

TA  
Set to "arbitrary" (default in waitingtime_regseq_whalevolution.m)

qvec  
0.5

---

# Notes

The calculations include a preparation step computing the ratio between the expected value and the median of an exponential distribution. This is used when interpreting waiting-time results.
