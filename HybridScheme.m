
%%%%%%%%%%%%%%%%%%%%% The Hybrid scheme 2d %%%%%%%%%%%%%%%%%%%%%%%%%%

% plots a single realisation of a VMMA

% Claudio Heinrich, August 2016

clear all;
close all;



kappa=2;                %depth of the Hybrid scheme
a=-0.5;                 %roughness parameter
n=100;                  %grid resolution
g=0.2;                  %the parameter gamma 
N = floor(n^(1+g));     %the integral range is N/n



%% The volatility factor sigma

%load('volatility') % needs to be square matrix of length 2(floor(n^(1+g))+n)+1 
sigma=vol(n,N);

%% The matrix containing the evaluation locations \|b_k\|

bMat= bMatSimple(N);


%% Auxiliary objects
%LgMat contains the values L(\|k\|/n) for k\in {-kappa,...,kappa}^2, and
%the values g(k/n) for k\in {-N,...,N}^2 - {-kappa,...,kappa}^2

LgMat=LgMatMatern(n,N,a,kappa,bMat);
%LgMat=LgMatExponential(n,N,a,kappa,bMat);


%% Simulate Gaussian RVs

C=Cov(kappa,a,n);

W0=mvnrnd(zeros((2*kappa+1)^2+1,1),C,(2*n+2*kappa+1)^2).';
W01=reshape(W0(1,:),[2*n+2*kappa+1,2*n+2*kappa+1]);
% Contains the W^n_{\bi} for \bi in {-n-kappa,...,n+kappa}^2
W02=reshape(W0(2:end,:),[2*kappa+1,2*kappa+1,2*n+2*kappa+1,2*n+2*kappa+1]); 
% Contains the W^n_{\bi,\bj} for \bi in {-n-kappa,...,n+kappa}^2

We1=normrnd(0,1/n^2,[2*n+2*N+1,2*n+2*N+1]);
We1(N-kappa+2:N+2*n+kappa+2,N-kappa+2:N+2*n+kappa+2)=W01;
% We1(i1+N+n+1,i2+N+n+1) contains W^n_{\bi} with \bj=(j1,j2)
% and \bi = (i1,i2) for i1,i2 in {-n-N,n+N}^2

%%
% for creating samples based on the same white noise, e.g. to demonstrate the 
% impact of different volatilities, comment the definitions of W0,W01,W02 and We1
% above and load the corresponding file here:

% load('WNexp');

%%

We2=zeros(2*kappa+1,2*kappa+1,2*n+2*N+1,2*n+2*N+1);
We2(:,:,N-kappa+2:N+2*n+kappa+2,N-kappa+2:N+2*n+kappa+2)=W02;
%4d array, We2(j1+kappa+1,j2+kappa+1,i1+N+n+1,i2+N+n+1) contains W^n_{\bi,\bj} with \bj=(j1,j2)
% and \bi = (i1,i2) for i1,i2 in {-n-kappa,n+kappa}^2, and zeros else     



%% Simulation of \tilde X, i.e. the integral around 0 

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
        B=LgMat(N-kappa+1:N+kappa+1,N-kappa+1:N+kappa+1).*Wshift(:,:,i1+n+1,i2+n+1);
        X1(i1+n+1,i2+n+1)= sum(B(:));
    end
end



%% Simulation of \hat X, that is the integral away from 0

% gMat contains the values g(\bk/n) for \bk\in\{-N,...,N\}^2\setminus\{-kappa,...,kappa\}^2  
% it contains 0 at the positions corresponding to \{-kappa,...,kappa\}^2

gMat=LgMat;
gMat(N+1-kappa:N+1+kappa,N+1-kappa:N+1+kappa)=zeros(2*kappa+1);

 
X2=conv2fft(sigma.*We1,gMat,'valid');
 
X=X1+X2;

%% plotting



surf(-1:1/n:1,-1:1/n:1,X,'EdgeColor','none');

set(gca,'FontSize',12)
xlabel('$t_1$','Interpreter','latex')
ylabel('$t_2$','Interpreter','latex')
zlabel('$X_{\bf{t}}$','Interpreter','latex')
title(['$\alpha=$' num2str(a) ],'interpreter','latex','FontSize',28)

