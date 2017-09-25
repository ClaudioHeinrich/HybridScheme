function [ b ] = bMatSimple( N )

% Contains the evaluation points \bb_{\bj}

b=zeros(2*N+1);



for i=-N:N
    for j=-N:N
        b(i+N+1,j+N+1)=norm([i,j]);    
    end
end

% 
% for i=0:N
%     for j=0:i
%         b(i+N+1,j+N+1)=norm([i,j]);    
%     end
% end
% 
% b(N+1,N+1)=0; %b_(0,0)
% 
% 
% for i=0:N-1
%     for j=i+1:N
%         b(i+N+1,j+N+1)=b(j+N+1,i+N+1);
%     end
% end
% 
% for i=-N:-1
%     for j=-N:-1
%          b(i+N+1,j+N+1)=b(-i+N+1,-j+N+1);
%     end
%     for j=0:N
%          b(i+N+1,j+N+1)=b(-i+N+1,j+N+1);
%     end
% end
% 
% for i=0:N
%     for j=-N:-1
%          b(i+N+1,j+N+1)=b(i+N+1,-j+N+1);
%     end
% end
% 


end

