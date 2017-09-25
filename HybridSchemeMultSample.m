%%%%%%%%%%%%%%%%%%%%% The Hybrid scheme 2d %%%%%%%%%%%%%%%%%%%%%%%%%%


% generates Monte-Carlo samples for a range of parameters alpha, kappa and
% n. Output is a 6-dimensional array with the first dimension indexing the
% Monte Carlo sample, the second n, the third kappa, the fourth alpha, and
% the last two dimensions containing the 2n+1 x 2n+1 values of the process.

% Claudio Heinrich, August 2016

clear all;
close all;

sam=100;            % numbers of i.i.d. MC samples 
avec=[-.8:.2:-.2];  % vector containing the roughness parameters
kappavec=[1,2];     % vector containing the values for kappa
nvec=[75];          % vector containing the grid resolutions

g=0.3;    Nvec=floor(nvec.^(g+1));   
                    % parameter gamma and boundary for the integral
                            
X=nan(sam,length(nvec),length(kappavec),length(avec),2*max(nvec)+1,2*max(nvec)+1);

for nind=1:length(nvec)
    
 n=nvec(nind);
 
 N = Nvec(nind);  

%% Volatility

sigma=vol(n,N);

%% Matrix containing \|\bb_\bj\|

 bMat=bMatSimple(N);
 

for aind=1:length(avec)
    a=avec(aind);

    
for kappaind=1:length(kappavec)
 kappa=kappavec(kappaind);
 
 Mat=LgMatMatern(n,N,a,kappa,bMat);
 %Mat=LgMatExponential(n,N,a,kappa,bMat);

%% Simulate Gaussian RVs

C=Cov(kappa,a,n);

for samplecounter = 1:sam

%Generate Gaussian variables
W0=mvnrnd(zeros((2*kappa+1)^2+1,1),C,(2*n+2*kappa+1)^2).';

%Get the W^n_{i} for i in {-n-kappa,...,n+kappa}^2 as matrix
W01=reshape(W0(1,:),[2*n+2*kappa+1,2*n+2*kappa+1]);

%Get the W^n_{i,j} for i,j in {-n-kappa,...,n+kappa}^2 as array
W02=reshape(W0(2:end,:),[2*kappa+1,2*kappa+1,2*n+2*kappa+1,2*n+2*kappa+1]); 


We1=normrnd(0,1/n^2,[2*n+2*N+1,2*n+2*N+1]);
We1(N-kappa+2:N+2*n+kappa+2,N-kappa+2:N+2*n+kappa+2)=W01;
% We1(i1+N+n+1,i2+N+n+1) contains W^n_{i} with j=(j1,j2)
% and i = (i1,i2) for i1,i2 in {-n-N,n+N}^2


We2=zeros(2*kappa+1,2*kappa+1,2*n+2*N+1,2*n+2*N+1);
We2(:,:,N-kappa+2:N+2*n+kappa+2,N-kappa+2:N+2*n+kappa+2)=W02;
%4d array, We2(j1+kappa+1,j2+kappa+1,i1+N+n+1,i2+N+n+1) contains W^n_{i,j} with j=(j1,j2)
% and i = (i1,i2) for i1,i2 in {-n-kappa,n+kappa}^2, and zeros else     



%% Simulation of \tilde X

%Wshift contains sigma_{\bi-\bk}W_{\bi-\bk,\bk} at position (k1+kappa+1,k2+kappa+1,i1+n+1,i2+n+1) 

Wshift=nan(2*kappa+1,2*kappa+1,2*n+1,2*n+1);
for k1=-kappa:kappa
    for k2=-kappa:kappa
        for i1=-n:n
            for i2=-n:n
                Wshift(k1+kappa+1,k2+kappa+1,i1+n+1,i2+n+1)=sigma(i1-k1+N+n+1,i2-k2+N+n+1)*We2(k1+kappa+1,k2+kappa+1,i1-k1+N+n+1,i2-k2+N+n+1);
            end
        end
    end
end


X1=nan(2*n+1);

for i1= -n:n
    for i2 = -n:n
        B=Mat(N-kappa+1:N+kappa+1,N-kappa+1:N+kappa+1).*Wshift(:,:,i1+n+1,i2+n+1);
        X1(i1+n+1,i2+n+1)= sum(B(:));
    end
end


%% Simulation of \hat X

% gMat contains the values g(\bk/n) for \bk\in\{-N,...,N\}^2\setminus\{-kappa,...,kappa\}^2  
% it contains 0 at the positions corresponding to \{-kappa,...,kappa\}^2

gMat=Mat;
gMat(N+1-kappa:N+1+kappa,N+1-kappa:N+1+kappa)=zeros(2*kappa+1);

 sigWe1=sigma.*We1;
 
 X2=conv2fft(sigWe1,gMat,'valid');

 X(samplecounter,nind,kappaind,aind,1:2*n+1,1:2*n+1)=X1+X2;

end
end
end
end

