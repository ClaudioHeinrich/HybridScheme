function [  Mat ] = LgMatExponential(n,N,a,kappa,b)

%The output matrix contains the values L(\bk/n) for \bk\in {-kappa,...,kappa}^2, and
%the values g(\bk/n) for \bk\in {-N,...,N}^2 - {-kappa,...,kappa}^2 for the
%case L(x)=exp(-x)
%In order to minimize function calls, compute LgMat only on half a 
%quadrant and exploit isometry.


%exponential L:
 Lfct=@(x)(exp(-norm(x)));
 Lfct1d =@(x)(exp(-abs(x))) ;

 %value of L at the optimal discretisation point in [-1/n,1/n]^2. The
 %numeric integration becomes degenerate when a<-0.8.
if a < -0.8
    optL=1;
else
    intfct = @(x)(Lfct1d(x/n).*(x.^(2*a+1)).*(pi/4-(x>=1/2).*acos(sqrt(2)*x)));
    optL =integral(intfct,0,1/sqrt(2))./TriInt0(1/2,a);
end
    



Mat=nan(2*N+1);
for i=0:N
    for j=0:i
        if abs(i)>kappa 
            Mat(i+N+1,j+N+1)=Lfct1d(b(i+N+1,j+N+1)/n)*(b(i+N+1,j+N+1)/n)^a;
        else if abs(i)>0
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

