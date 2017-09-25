%returns the covariance matrix of the vector (W_0^n,W^1',...,W^(2kappa+1)'), 
%where W^k'=(W_{(0,0)',(-kappa,k-kappa-1)'},...,W_{(0,0)',(kappa,k-kappa-1)'})

function  C = Cov(kappa,a,n)

A=Cov2(kappa,a,n);
A=reshape(A,[(2*kappa+1)^2,(2*kappa+1)^2]);

B=Cov1(kappa,a,n);
B=reshape(B,[], 1);

C=nan((2*kappa+1)^2+1);
C(1,1)=1/(n^2);
C(2:(2*kappa+1)^2+1,1)=B;
C(1,2:(2*kappa+1)^2+1)=B.';
C(2:(2*kappa+1)^2+1,2:(2*kappa+1)^2+1)=A;
end

%%%%%%%%%%%%%

function C1 = Cov1(kappa,a,n)


%returns the covariances C_{1,\bj}. The output matrix C1 is defined as C1(j,k)=C_{1,(j-kappa-1,k-kappa-1)'}

TriMa=TriIntMat(kappa,a/2);

C1=nan(2*kappa+1,2*kappa+1);

C1(kappa+1,kappa+1)=8*TriMa(1,1);


% C_{(j,j)',(j,j)'}, j > 0


coor=nan(1,4);

for j=1:kappa
    value=2*TriMa(j+1,j+1);
    
    coor(1)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,j+kappa+1);
    coor(2)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,j+kappa+1);
    coor(3)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,-j+kappa+1);
    coor(4)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,-j+kappa+1);
    
    C1(coor)=value;

end


% C_{1,(1,0)'}


if kappa>0
    value=2*(TriMa(2,1)-TriMa(2,2)-TriMa(1,1));
    
    coor(1)=sub2ind([2*kappa+1 2*kappa+1],kappa+2,kappa+1);
    coor(2)=sub2ind([2*kappa+1 2*kappa+1],kappa+1,kappa+2);
    coor(3)=sub2ind([2*kappa+1 2*kappa+1],kappa,kappa+1);
    coor(4)=sub2ind([2*kappa+1 2*kappa+1],kappa+1,kappa);
    
    C1(coor)=value;
end


% C_{1,(j,0)'}, j > 1

if kappa > 1
    for j=2:kappa
        value=2*(TriMa(j+1,1)-TriMa(j+1,2)-TriMa(j,1)+TriMa(j,2));
        
        coor(1)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,kappa+1);
        coor(2)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,kappa+1);
        coor(3)=sub2ind([2*kappa+1 2*kappa+1],kappa+1,j+kappa+1);
        coor(4)=sub2ind([2*kappa+1 2*kappa+1],kappa+1,-j+kappa+1);
    
        C1(coor)=value;
    end
end


% C_{1,(j,k)'},  0 < k = j+-1

coor=nan(8,1);

if kappa>1
    for j=2:kappa
         value=TriMa(j+1,j)-TriMa(j+1,j+1)-TriMa(j,j);
                  
         coor(1)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,j+kappa);
         coor(2)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,j+kappa);
         coor(3)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,-j+kappa+2);
         coor(4)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,-j+kappa+2);
         coor(5)=sub2ind([2*kappa+1 2*kappa+1],j+kappa,j+kappa+1);
         coor(6)=sub2ind([2*kappa+1 2*kappa+1],j+kappa,-j+kappa+1);
         coor(7)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+2,j+kappa+1);
         coor(8)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+2,-j+kappa+1);
         
      
         C1(coor)=value;
    end
end


% C_{1,(j,k)'},  0 < k < j-1

if kappa>2
    for j=3:kappa
        for k=1:j-2
            value=TriMa(j+1,k+1)-TriMa(j+1,k+2)-TriMa(j,k+1)+TriMa(j,k+2);
            coor(1)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,k+kappa+1);
            coor(2)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,k+kappa+1);
            coor(3)=sub2ind([2*kappa+1 2*kappa+1],j+kappa+1,-k+kappa+1);
            coor(4)=sub2ind([2*kappa+1 2*kappa+1],-j+kappa+1,-k+kappa+1);
            coor(5)=sub2ind([2*kappa+1 2*kappa+1],k+kappa+1,j+kappa+1);
            coor(6)=sub2ind([2*kappa+1 2*kappa+1],k+kappa+1,-j+kappa+1);
            coor(7)=sub2ind([2*kappa+1 2*kappa+1],-k+kappa+1,j+kappa+1);
            coor(8)=sub2ind([2*kappa+1 2*kappa+1],-k+kappa+1,-j+kappa+1);
            
            
            C1(coor)=value;
        end
    end
end
C1=n^(-2-a)*C1;

end

%%%%%%%%%%%%%%%

function covM = Cov2(kappa,a,n)

%Returns the (2kappa+1)^4 array covM with entries
%covM(j1,j2,k1,k2)= C_{(j_1-kappa-1,j_2-kappa-1)',(k_1-kappa-1,k_2-kappa-1)'}.


format long;
% The matrix TriIntMat contains integrals of \|x\|^a over triangular sets
% as described in the covariance section in the paper

TriMa=TriIntMat(kappa,a);

%% Create (2kappa+1)^4 array for the covariance entries. 

covM = nan(2*kappa+1,2*kappa+1,2*kappa+1,2*kappa+1);

% C_{(0,0),(0,0)}

covM(kappa+1,kappa+1,kappa+1,kappa+1)=8*TriMa(1,1);

% C_{(j,j),(j,j)}, j > 0

for j=1:kappa
    coor=symind([j+kappa+1; j+kappa+1; j+kappa+1; j+kappa+1],kappa);
    covM(coor)=2*TriMa(j+1,j+1);
end


% C_{(1,0),(1,0)}
if kappa>0
    value=2*(TriMa(2,1)-TriMa(2,2)-TriMa(1,1));
    coor=symind([kappa+2; kappa+1; kappa+2; kappa+1],kappa);
    covM(coor)=value;
end


% C_{(j,0),(j,0)}, j > 1

if kappa > 1
    for j=2:kappa
        value=2*(TriMa(j+1,1)-TriMa(j+1,2)-TriMa(j,1)+TriMa(j,2));
        coor=symind([j+kappa+1; kappa+1; j+kappa+1; kappa+1],kappa);
        covM(coor)=value;
    end
end


% C_{(j,k),(j,k)},  0 < k = j-1

if kappa>1
    for j=2:kappa
         value=TriMa(j+1,j)-TriMa(j+1,j+1)-TriMa(j,j);
         coor=symind([j+kappa+1; j+kappa; j+kappa+1; j+kappa],kappa);
         covM(coor)=value;
    end
end


% C_{(j,k),(j,k)},  0 < k < j-1

if kappa>2
    for j=3:kappa
        for k=1:j-2
            value=TriMa(j+1,k+1)-TriMa(j+1,k+2)-TriMa(j,k+1)+TriMa(j,k+2);
            coor=symind([j+kappa+1; k+kappa+1; j+kappa+1; k+kappa+1],kappa);
            covM(coor)=value;
        end
    end
end

% C_{(j1,j2),(k1,k2)}, (j1,j2) neq 0, (k1,k2) neq 0, (j1,j2) neq (k1,k2):
% Filling in the missing covariances by numeric integration

for j1=1:kappa
    for j2=0:j1
        for k1=-kappa:kappa
            for k2=-kappa:kappa
                if isnan(covM(j1+kappa+1,j2+kappa+1,k1+kappa+1,k2+kappa+1))
                    fun=@(x,y)(((j1-x).^2+(j2-y).^2).^(a/2).*((k1-x).^2+(k2-y).^2).^(a/2));
                    value=integral2(fun,-0.5,0.5,-0.5,0.5,'AbsTol',1e-10,'RelTol',0);
                    coor=symind([j1+kappa+1; j2+kappa+1; k1+kappa+1; k2+kappa+1],kappa);
                    covM(coor)=value;
                end
            end
        end
    end
end

covM=n^(-2-2*a)*covM;

end


%%%%%%%%%%%%%%%%

function x = TriangleIntegral(j1,j2,a) 


%Takes input (j1,j2) with 0<j2<j1 and computes the
% integral of \|x\|^2a over the set {(x,y) : 0 < x < i-1/2, j-1.5 < y < x},
% see Appendix B for details 

x=((j1.^(2*a+2)+j2.^(2*a+2))/(2^(1.5)*(a+1))).*hypergeom([0.5,1.5+a],1.5,0.5);
x=x-(j1.*j2.^(2*a+2)).*hypergeom([0.5,1.5+a],1.5,j1.^2/(j1.^2+j2.^2))./(2*(a+1)*sqrt(j1.^2+j2.^2));
x=x-(j1.^(2*a+2).*j2).*hypergeom([0.5,1.5+a],1.5,j2.^2/(j1.^2+j2.^2))/(2*(a+1)*sqrt(j1.^2+j2.^2));
end

%%%%%%%%%%%%%%%%%

function A = TriIntMat(kappa,a)

% TriIntMat contains the integrals of \|x\|^a over 
% triangular sets as follows.  
% TriIntMat is a symmetric matrix and witht the entries in the first column
% TriIntMat(i,1) = integral of \|x\|^a over the set {(x,y) : 0 < x < i-1/2, 0 < y < x} 
% and for all other columns, i.e. with j>1
% TriIntMat(i,j) = integral of \|x\|^a over the set {(x,y) : 0 < x < i-1/2, j-1.5 < y < x} 
% See Appendix B for details.

A=nan(kappa+1,kappa+1);

A(1,1)= TriInt0(0.5,a);

for j1=2:kappa+1
A(1,j1)=TriInt0(j1-0.5,a);
A(j1,1)=A(1,j1);
    for j2=2:j1
        A(j1,j2)=TriangleIntegral(j1-0.5,j2-1.5,a);
        A(j2,j1)=A(j1,j2);
    end
end
end

%%%%%%%%%%%%%%%%%

function Y = symind(A,kappa) 

% The argument A is a vector with 4 entries between 1 and 2kappa+1 
% representing a coordinate in the 4 dimensional covariance array of size 
% (2kappa+1)^4 containing C_{(i_1,j_1),(i_2,j_2)}.
% The function returns the 16 coordinates of the covariance array 
% that have the same entry as coordinate A by symmetry arguments.
% The output is a 1 x 16 rowvector containing these coordinates as linear
% indices

% we introduce the array M, containing on each page a 4 x 4 matrix
% corresponding to one of the 8 spatial symmetries that map our grid onto
% itself

M2= [-1 0 0 0;0 1 0 0;0 0 -1 0;0 0 0 1];
M3= [0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0];
M4= [0 -1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0];

M=cat(3,eye(4),-eye(4),M2,-M2,M3,-M3,M4,-M4);

% The matrix S corresponds to the symmetry C_{j,k}=C_{k,j}

S=[0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];

%coordinate shift
A=A-kappa-1;

X=nan(4,16);
for i=1:8
    X(:,i)= M(:,:,i)*A+kappa+1;
end
A=S*A;
for i=9:16
    X(:,i)= M(:,:,i-8)*A+kappa+1;
end

Y=nan(1,16);
for j=1:16
    Y(j)=sub2ind([2*kappa+1 2*kappa+1 2*kappa+1 2*kappa+1],X(1,j),X(2,j),X(3,j),X(4,j));
end
end


