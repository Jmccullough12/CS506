%CS506 Project Test Code
%Jeb Kilfoye, Jacob McCullough
%You can enable plotting, pass in your desired number of agents n.
%Various properties of the calculations can be seen by unsurpressing output 
%Interesting ones include max_distances, total_distance, robots_x, robots_y

function h = circle_n_agents(n)
x = 0;
y = 0;
r = 1;

d = r*2;
px = x-r;
py = y-r;

%turn plotting off to speed up calculations
plotting = true;

if plotting
    clf;
    hold on;
    h = rectangle('Position',[px py d d],'Curvature',[1,1]);
    xlim([-1,1]);
    ylim([-1,1]);
end


s = random_circle_point();
sx = s(1);
sy = s(2);
if plotting
    plot(sx,sy,'r.');
end

angle_offset_degrees = 360/n;
robots_x = zeros(n,1);
robots_y = zeros(n,1);
perturb = .001;

epsilon = perturb;
D = 2; %space diameter
G = 1; %maximum gradient
lambda = 1/exp(-2); %1-e^(-x) for x in [0,2] has minima grad at x=2, one over that
T = D^2*G^2/epsilon^2; %max steps
eta = D/(G*sqrt(T)); %step size


%random vs evenly distributed
for i = 1:n
   %rob_i = random_circle_point(); 
   %robots_x(i) = rob_i(1);
   %robots_y(i) = rob_i(2);
   robots_x(i) = (1/2)*cosd(0+(i-1)*angle_offset_degrees)+rand(1)*2*perturb-perturb;
   robots_y(i) = (1/2)*sind(0+(i-1)*angle_offset_degrees)+rand(1)*2*perturb-perturb;
end

if plotting
    for i = 1:n
        plot(robots_x(i),robots_y(i), 'b.','MarkerSize',10);
        text(robots_x(i),robots_y(i),string(i));
    end
end

concentrations = zeros(n,1);
for i = 1:n
    concentrations(i) = concentration(robots_x(i),robots_y(i),sx,sy);
end
concentrations;
j_delta = .001;
gradients_x = zeros(n,1);
gradients_y = zeros(n,1);
for i = 1:n
    result = get_gradient(robots_x(i),robots_y(i),sx,sy,j_delta);
    gradients_x(i) = result(1);
    gradients_y(i) = result(2);
end
gradients = [gradients_x, gradients_y];
perp_slopes = zeros(n,1);
perp_intercepts = zeros(n,1);
for i = 1:n
    x1 = robots_x(i);
    y1 = robots_y(i);
    m1=-(1/(gradients_y(i)/gradients_x(i)));
    perp_slopes(i) = m1;
    perp_intercepts(i) = y1-m1*x1;
    x = linspace(-1,1,50);
    y = perp_slopes(i)*(x-x1)+y1;
    if plotting 
        plot(x,y);
    end
end

closed_convex = true;
done = false;
for i = 1:n
    for j = 1:n
        if j ~= i
            i;
            j;
            s1 = sign(perp_slopes(i)*robots_x(j)+perp_intercepts(i)-robots_y(j));
            s2 = sign(perp_slopes(i)*sx+perp_intercepts(i)-sy);
            if  s1 ~= s2  
                closed_convex = false;
                done = true;
                break;
            end
        end
        if done
            break;
        end
    end
end
closed_convex;

%if closed, find intersections, take convex hull, find max dist point on
%hull

hull = [-1, -1;
        -1,  1;
         1,  1;
         1, -1;
        -1, -1];
for i = 1:n
    if sy > perp_slopes(i)*sx+perp_intercepts(i)
        hull = cutpolygon(hull,[robots_x(i),robots_y(i);0,perp_intercepts(i)],'B');
    else
        hull = cutpolygon(hull,[robots_x(i),robots_y(i);0,perp_intercepts(i)],'T');
    end
end
hull;
if plotting
   plot(hull(:,1),hull(:,2),'r') 
end



%if open, find lines used in source segment, find max dist between source
%and circular perimeter, and line segment intersections
max_distances = zeros(n,1);
for i = 1:n
    max_distances(i) = most_distant_point(hull,closed_convex,robots_x(i),robots_y(i),perp_slopes,perp_intercepts);
end

total_distance = 0;
active_robots = ones(n,1);
for i = 1:T-1

    %check if best case for robot j is worse than worst case for robot k
    %using gradient bound
%     for j = 1:n
%         for k = 1:n
%             max_dist_bound(T-i,epsilon,lambda, G, ...
%                     most_distant_point(hull,0,robots_x(j),robots_y(j), ...
%                     perp_slopes,perp_intercepts))
%             if j ~= k & max_dist_bound(T-i,epsilon,lambda, G, ...
%                     most_distant_point(hull,0,robots_x(j),robots_y(j), ...
%                     perp_slopes,perp_intercepts)) < distance_to_hull(hull,[robots_x(k),robots_y(k)])
%                 active_robots(k)=false;
%             end
%         end
%     end

    %check if best case for robot j is worse than worst case for robot k
    %using lin dist

    for j = 1:n
        for k = 1:n
            if j ~= k & most_distant_point(hull,0,robots_x(j),robots_y(j), ...
                    perp_slopes,perp_intercepts) < distance_to_hull(hull,[robots_x(k),robots_y(k)])
                active_robots(k)=false;
            end
        end
    end



    step = gradient_descent_step(robots_x,robots_y,gradients(:,1),gradients(:,2),eta,active_robots);
    for j = 1:n
        total_distance = total_distance + distance(step(j,1),step(j,2),robots_x(j),robots_y(j));
    end
    if mod(i,20)==1 & plotting
        for j = 1:n
            if active_robots(j)
                plot(step(j,1),step(j,2), 'b.','MarkerSize',10);
                text(step(j,1),step(j,2),strcat(string(j),'*'));
            end
        end
    end
    robots_x=step(:,1);
    robots_y=step(:,2);
    gradients_x = zeros(n,1);
    gradients_y = zeros(n,1);
    for i = 1:n
        result = get_gradient(robots_x(i),robots_y(i),sx,sy,j_delta);
        gradients_x(i) = result(1);
        gradients_y(i) = result(2);
    end
    gradients = [gradients_x, gradients_y];
    perp_slopes = zeros(n,1);
    perp_intercepts = zeros(n,1);
    for i = 1:n
        x1 = robots_x(i);
        y1 = robots_y(i);
        m1=-(1/(gradients_y(i)/gradients_x(i)));
        perp_slopes(i) = m1;
        perp_intercepts(i) = y1-m1*x1;
    end
    closed_convex = true;
    done = false;
    for i = 1:n
        for j = 1:n
            if j ~= i
                i;
                j;
                s1 = sign(perp_slopes(i)*robots_x(j)+perp_intercepts(i)-robots_y(j));
                s2 = sign(perp_slopes(i)*sx+perp_intercepts(i)-sy);
                if  s1 ~= s2  
                    closed_convex = false;
                    done = true;
                    break;
                end
            end
            if done
                break;
            end
        end
    end
    closed_convex;
    
%     if closed, find intersections, take convex hull, find max dist point on
%     hull
    
    hull = [-1, -1;
            -1,  1;
             1,  1;
             1, -1;
            -1, -1];
    for i = 1:n
        if active_robots(i)
            if sy > perp_slopes(i)*sx+perp_intercepts(i)
                hull = cutpolygon(hull,[robots_x(i),robots_y(i);0,perp_intercepts(i)],'B');
            else
                hull = cutpolygon(hull,[robots_x(i),robots_y(i);0,perp_intercepts(i)],'T');
            end
        end
    end
    max_distances = zeros(n,1);
    for i = 1:n
        max_distances(i) = most_distant_point(hull,closed_convex,robots_x(i),robots_y(i),perp_slopes,perp_intercepts);
    end

    %done checks if any robot is close enough to the source
    done = zeros(n,1);
    for j = 1:n
        if distance(robots_x(j),robots_y(j),sx,sy)<=2*epsilon
            done(j) = 1;
        end
    end
    if sum(done)==1
        break;
     end
 end
 total_distance


if plotting
    hold off;
    daspect([1,1,1])
end

end

%linear distance finder
function dist = distance(x1,y1,x2,y2)
dist = sqrt(((x1-x2)^2)+((y1-y2)^2));
end

%Our f(x) function for the concentration
function conc = concentration(x,y,sx,sy)
conc = 1-exp(-distance(x,y,sx,sy));
end

%samples gradient by taking a delta step and comparing differences
function grad = get_gradient(x,y,sx,sy,delta)
temp_x = x + delta;
temp_y = y + delta;
grad_x = (concentration(temp_x,y,sx,sy)-concentration(x,y,sx,sy))/delta;
grad_y = (concentration(x, temp_y,sx,sy)-concentration(x,y,sx,sy))/delta;
grad = [grad_x,grad_y];
end

%random point in the unit circle
function r_point = random_circle_point
x = -5;
y = -5;
while sqrt(x^2+y^2)>1
    x = 2*rand(1)-1;
    y = 2*rand(1)-1;
end
r_point = [x,y];
end

%calculates next step of gradient descent
function new_pos = gradient_descent_step(rx,ry, grads_x,grads_y, eta,active_robots)
    new_pos = zeros(length(rx),2);
    n_grad = zeros(length(rx),2);
    for i=1:length(rx)
        if active_robots(i)
            n_grad(i,:) = -1*([grads_x(i),grads_y(i)])/norm([grads_x(i),grads_y(i)]);
            new_pos(i,:) = [rx(i)+eta*n_grad(i,1),ry(i)+eta*n_grad(i,2)];
        end
    end
end

%bound on max distance considering a non-straightline path to minima
function dist_bound = max_dist_bound(iters_left, epsilon, grad_lower,grad_upper,lin_dist)
    a = grad_lower*(exp(-lin_dist)-exp(-epsilon));
    b = (grad_upper^2/epsilon)*((grad_lower^2)*(1-exp(-lin_dist))^2-epsilon^2);
    c = (iters_left-1)*epsilon;
    dist_bound = a+b+c;
end

%maax linear distance from a point to any point in the convex hull
function max_dist = most_distant_point(hull, closed_hull, rx,ry,slopes,intercepts)
    p = [inf,inf]; %max dist point for printout
        max_dist = 0; 
        for i = 1:length(hull) %check points in convex hull if within radius 1
            if distance(hull(i,1),hull(i,2),0,0) <= 1
                cur_dist = distance(hull(i,1),hull(i,2),rx,ry);
                if cur_dist > max_dist
                    max_dist = cur_dist;
                    p=[hull(i,1),hull(i,2)];
                end
            end
        end
        x1 = sqrt(rx^2/(rx^2+ry^2));%check antipodal point
        y1 = ry*x1/rx;
        x2 = -x1;
        y2 = ry*x2/rx;
        d1 = distance(x1,y1,rx,ry);
        d2 = distance(x2,y2,rx,ry);
        if inhull([x1,y1],hull,[] ,.01) & d1 > max_dist
            max_dist = d1;
            p=[x1,y1];
        elseif inhull([x2,y2],hull,[],.01) & d2 > max_dist
            max_dist = d2;
            p=[x2,y2];
        end
        for i = 1:length(slopes) %Find intercepts with circle x^2+y^2=1
            m = slopes(i);
            b = intercepts(i);
            x1 = (-2*m*b+sqrt(4*(m^2)*(b^2)-4*(m^2+1)*(b^2-1)))/(2*(m^2+1));
            y1 = m*x1+b;
            x2 = (-2*m*b-sqrt(4*(m^2)*(b^2)-4*(m^2+1)*(b^2-1)))/(2*(m^2+1));
            y2 = m*x2+b;
            d1 = distance(x1,y1,rx,ry);
            d2 = distance(x2,y2,rx,ry);
            if inhull([x1,y1],hull,[],.01) & d1 > max_dist
                max_dist = d1;
                p=[x1,y1];
            elseif inhull([x2,y2],hull,[],.01) & d2 > max_dist
                max_dist = d2;
                p=[x2,y2];
            end
        end
    %end
    
end

%Find the shortest distance to the convex hull for a point p
function dist = distance_to_hull(hull,p)
    n=size(hull,1);
    dist = Inf;
    for i = 1:n-1
        a = horzcat(hull(i,:)-hull(i+1,:),[0]);
        b = horzcat(p - hull(i+1,:),[0]);
        d = norm(cross(a,b))/norm(a);
        d1 = distance(hull(i,1),hull(i,2),p(1),p(2));
        d2 = distance(hull(i+1,1),hull(i+1,2),p(1),p(2));
        if max(d1,d2)< d & d < dist
                dist = d;
        elseif min(d1,d2) < dist
            dist = min(d1,d2);
        end
    end
end

