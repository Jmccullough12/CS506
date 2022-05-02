function h = circle2(x,y,r,n)
clf;
d = r*2;
px = x-r;
py = y-r;
hold on;
h = rectangle('Position',[px py d d],'Curvature',[1,1]);
s = random_circle_point();
sx = s(1);
sy = s(2);
plot(sx,sy,'r.');
angle_offset_degrees = 360/n;
robots_x = zeros(n,1);
robots_y = zeros(n,1);
perturb = .001;
for i = 1:n
   rob_i = random_circle_point()
   robots_x(i) = rob_i(1);
   robots_y(i) = rob_i(2);
   %robots_x(i) = (1/2)*cosd(0+(i-1)*angle_offset_degrees)+rand(1)*2*perturb-perturb;
   %robots_y(i) = (1/2)*sind(0+(i-1)*angle_offset_degrees)+rand(1)*2*perturb-perturb;
end

for i = 1:n
    plot(robots_x(i),robots_y(i), 'b.','MarkerSize',10);
    text(robots_x(i),robots_y(i),string(i));
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
    plot(x,y);
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
hull
plot(hull(:,1),hull(:,2),'r')


%if open, find lines used in source segment, find max dist between source
%and circular perimeter, and line segment intersections
max_distances = zeros(n,1);
for i = 1:n
    max_distances(i) = most_distant_point(hull,closed_convex,robots_x(i),robots_y(i),perp_slopes,perp_intercepts);
end

xlim([-1,1]);
ylim([-1,1]);
hold off;
daspect([1,1,1])

end

function dist = distance(x1,y1,x2,y2)
dist = sqrt(((x1-x2)^2)+((y1-y2)^2));
end

function conc = concentration(x,y,sx,sy)
conc = 1-exp(-distance(x,y,sx,sy));
end

function grad = get_gradient(x,y,sx,sy,delta)
temp_x = x + delta;
temp_y = y + delta;
grad_x = (concentration(temp_x,y,sx,sy)-concentration(x,y,sx,sy))/delta;
grad_y = (concentration(x, temp_y,sx,sy)-concentration(x,y,sx,sy))/delta;
grad = [grad_x,grad_y];
end

function r_point = random_circle_point
x = -5;
y = -5;
while sqrt(x^2+y^2)>1
    x = 2*rand(1)-1;
    y = 2*rand(1)-1;
end
r_point = [x,y];
end

function max_dist = most_distant_point(hull, closed_hull, rx,ry,slopes,intercepts)
    if closed_hull
        max_dist = 0;
        for i = 1:length(hull)
            cur_dist = distance(hull(i,1),hull(i,2),rx,ry);
            if cur_dist > max_dist
                max_dist = cur_dist;
            end
        end
    else
        max_dist = 0;
        for i = 1:length(hull)
            if distance(hull(i,1),hull(i,2),0,0) <= 1
                cur_dist = distance(hull(i,1),hull(i,2),rx,ry);
                if cur_dist > max_dist
                    max_dist = cur_dist;
                end
            end
        end
        x1 = sqrt(rx^2/(rx^2+ry^2));
        y1 = ry*x1/rx;
        x2 = -x1;
        y2 = ry*x2/rx;
        d1 = distance(x1,y1,rx,ry);
        d2 = distance(x2,y2,rx,ry);
        if inhull([x1,y1],hull) & d1 > max_dist
            max_dist = d1;
        elseif inhull([x2,y2],hull) & d2 > max_dist
            max_dist = d2;
        end
        for i = 1:length(slopes)
            m = slopes(i);
            b = intercepts(i);
            x1 = (-2*m*b+sqrt(4*(m^2)*(b^2)-4*(m^2+1)*(b^2-1)))/(2*(m^2+1));
            y1 = m*x1+b;
            x2 = (-2*m*b-sqrt(4*(m^2)*(b^2)-4*(m^2+1)*(b^2-1)))/(2*(m^2+1));
            y2 = m*x2+b;
            if inhull([x1,y1],hull) & d1 > max_dist
                max_dist = d1;
            elseif inhull([x2,y2],hull) & d2 > max_dist
                max_dist = d2;
            end
        end
    end
end

