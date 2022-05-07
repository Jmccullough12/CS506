% CS-506 Homework3 - Problem 6
% Jeb and Jacob


data = [-1/2 1/3 1; 3/4 -1/4 1; 1/2 -1/3 -1; -3/4 1/4 -1];

n = 2;
T = 14;
eta = .22;
M = 3/4;
N = 2;
m = 4;

R1 = 1;

R2 = 2 * M;

u = ones(N+1,T+1);
U = ones(T,1);

loss = zeros(N+1,T);


w_t = ones(N,T);
for t=1:T
    
   
  
    U(t) = sum(u(:,t));
    for i=1:N
        w_t(i,t) = (R1 * u(i,t)) / U(t);
    end

    found = false;
    for j=1:m
        
        temp = dot(w_t(:,t), (data(j,1:N) * data(j,N+1)));
        if temp < 0
            for i=1:N
                loss(i,t) = -1 * (data(j,i) * data(j,N+1));
            end
            found = true;       
        end

        if found
            break;
        end
    end

    if ~found
        break;
    end
    
    for i=1:N+1
        u(i,t+1) = u(i,t) * (1 - eta * (loss(i,t) + R2/2)/R2);
        
    end
    
    
end

w_t(:,t)
sum(w_t(:,t))

