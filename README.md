# HybridScheme

implementation of the hybrid scheme for simulation of volatility modulated moving averages, 
as derived in the article available at arXiv:1709.01310.

The script file HybridScheme.m creates a plot of a VMMA, whereas HybridSchemeMultSample.m generates multiple Monte-Carlo 
samples of VMMAs for a range of different parameters.

The volatility process is specified in the file vol.m

As an example for a nontrivial volatility field we include the file volatility.mat, containing the volatility field that 
is shown in Figure 4 of the paper. In order to apply this example the parameters of the hybrid scheme need to be set to 
n=100 and g=0.2





