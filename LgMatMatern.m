function [  Mat ] = LgMatMatern(n,N,a,kappa,b)

%The output matrix contains the values L(k/n) for k\in {-kappa,...,kappa}^2, and
%the values g(k/n) for k\in {-N,...,N}^2 - {-kappa,...,kappa}^2 for the
%Matern covariance case.
%In order to minimize function calls, we compute LgMat only on half a quadrant and
%exploit spatial symmetries. 


%Matern covariance kernel:
 lambda=1; %scale parameter for Matern covariance, roughness parameter nu is a+1
 Lfct = @(x)( norm(x)^(-a/2)*besselk(a/2,lambda*norm(x))) ;
 Lfct1d = @(x)( abs(x).^(-a/2).*besselk(a/2,lambda*abs(x))) ;

%value of L at the optimal discretisation point in [-1/n,1/n]^2, see section 3 for details.
intfct = @(x)(Lfct1d(x/n).*(x.^(2*a+1)).*(pi/4-(x>=1/(2)).*acos(sqrt(2)*x)));
optL =integral(intfct,0,1/sqrt(2))./TriInt0(1/2,a);

 
 
Mat=nan(2*N+1);
for i=0:N
    for j=0:i
        if abs(i)>kappa | abs(j)>kappa
            Mat(i+N+1,j+N+1)=Lfct1d(b(i+N+1,j+N+1)/n)*(b(i+N+1,j+N+1)/n)^a;
        else
            Mat(i+N+1,j+N+1)=Lfct1d(b(i+N+1,j+N+1)/n);
        end
    end
end

%avoid singularity at 0
Mat(N+1,N+1) = optL;

% fill the rest of the matrix by using symmetries

for i=0:N-1
    for j=i+1:N
        Mat(i+N+1,j+N+1)=Mat(j+N+1,i+N+1);
    end
end

for i=-N:-1
    for j=-N:-1
         Mat(i+N+1,j+N+1)=Mat(-i+N+1,-j+N+1);
    end
    for j=0:N
         Mat(i+N+1,j+N+1)=Mat(-i+N+1,j+N+1);
    end
end

for i=0:N
    for j=-N:-1
         Mat(i+N+1,j+N+1)=Mat(i+N+1,-j+N+1);
    end
end


end

