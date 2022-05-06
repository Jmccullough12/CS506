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

x = zeros (1,N);

for t=1:T

   
    

    U(t) = sum(u(:,t));
    w_t = zeros(1,N);
    for i=1:N
        w_t(i) = (R1 * u(i,t)) / U(t);
    end
    
    found = false;
    for j=1:m
        
        temp = dot(w_t, data(j,1:N));
        if temp < 0
            for i=1:N+1
                loss(i,t) = -1 * data(j,i);
            end
            found = true;       
        end

        if found
            break;
        end
    end
    if ~found
        x = u(:,t+1);
        break;
    end


    for i=1:N+1
        u(i,t+1) = u(i,t) * (1 - eta * (loss(i,t) + R2/2)/R2);
    end
    
   
    
end

if sum(x) == 0
    x = u(:,T+1);
end

x

