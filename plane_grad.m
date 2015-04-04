function [ p_x,p_y,sz ] = plane_grad( p,rad,numx,numy )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sz = size(p);
p_x = zeros(sz(1),sz(2));
p_y = zeros(sz(1),sz(2));

radx = rad*numx;
rady = rad*numy; 

x=1:sz(1);
y=1:sz(2);

% interior gradients
for i=(1+radx):(sz(1)-radx)
    for j = (1+rady):(sz(2)-rady)
        p_x(i,j) = (p(i+radx,j) - p(i-radx,j))/(2.0*radx);
        p_y(i,j) = (p(i,j+rady) - p(i,j-rady))/(2.0*rady);
    end
end

% exterior x gradients
for i=1:radx
    for j = (1+rady):(sz(2)-rady)
        p_x(i,j) = (p(i+radx,j) - p(i,j))/(1.0*radx);
        p_y(i,j) = (p(i,j+rady) - p(i,j-rady))/(2.0*rady);
    end
end

for i=(sz(1)-radx+1):sz(1)
    for j = (1+rady):(sz(2)-rady)
        p_x(i,j) = (p(i,j) - p(i-radx,j))/(1.0*radx);
        p_y(i,j) = (p(i,j+rady) - p(i,j-rady))/(2.0*rady);
    end
end

%exterior y gradients
for i=(1+radx):(sz(1)-radx)
    for j = 1:rady
        p_x(i,j) = (p(i+radx,j) - p(i-radx,j))/(2.0*radx);
        p_y(i,j) = (p(i,j+rady) - p(i,j))/(1.0*rady);
    end
end

for i=(1+radx):(sz(1)-radx)
    for j = (sz(2)-rady+1):sz(2)
        p_x(i,j) = (p(i+radx,j) - p(i-radx,j))/(2.0*radx);
        p_y(i,j) = (p(i,j) - p(i,j-rady))/(1.0*rady);
    end
end

%xmax=max(max(p_x))
%ymax=max(max(p_y))

for i = 1:sz(1)
    for j = 1:sz(2)
    r = sqrt(p_x(i,j)^2 + p_y(i,j)^2);
    p_x(i,j) = p_x(i,j)/(.5+r);
    p_y(i,j) = p_y(i,j)/(.5+r);
    end
end


% save  'p_y.txt' p_y -ASCII;
% save  'p_x.txt' p_x -ASCII;
% 
% figure;
% 
% plot(y,p_y(11,:));
% pause(4.0)

% figure;
% surf(y,x,p_y);
% pause(10.0)
% corners are zero by default


% figure;
% spacing_x = 8;
% spacing_y = 6;
% starty = 5;
% quiver(y(starty:spacing_y:sz(2)),x(1:spacing_x:sz(1)),-p_y(1:spacing_x:sz(1),starty:spacing_y:sz(2)),-p_x(1:spacing_x:sz(1),starty:spacing_y:sz(2)));
% axis([0 sz(2) 0 sz(1)]);
% pause(2)
 
end

