function [ p,numx,numy,xdim,ydim ] = chair_Boeing777(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit)
% stuff to build a plane

% gridspacing
numx = 1;
numy = 1;
dx = 1.0/numx;
dy = 1.0/numy;


extension = 9; %for agents to exit faster

%plane dimensions
xdim = x_exit+space_front+2*rb*rwb+space_back+x_exit+space_front+2*re1*rwe+space_back+x_exit+space_front+2*re2*rwe+space_back+x_exit; %rows of airplane
ydim = 3*nce*cwe+2*ae; 
chairdim_b = ncb*cwb+extension;
chairdim_e = nce*cwe+extension;

%chair dimensions
seat_depth = sd*numx;
chair_depth = cd*numx;
legroom_b = lrb*numx;
chair_width_b = cwb*numy;
legroom_e = lre*numx;
chair_width_e = cwe*numy;

% Set the grid spacing
x = 0:dx:xdim;
y = 0:dy:ydim;
% grid spacing for chair part of plane
ycb = 0:dy:chairdim_b;
yce = 0:dy:chairdim_e;


% matrix size
M = length(x);
N = length(y);
% chair part of matrix
NCb = length(ycb);
NCe = length(yce);



p = zeros(M, N);

% i = row numbers #, x direction
% j = seat numbers #, y direction

% aisle function
for i=1:space_front+x_exit+rb*rwb
    for j=1:N
      p(i,j) = as*(x(i));%/xdim);   
    end
end
for i=(space_front+x_exit+rb*rwb+1):space_front+1.5*x_exit+2*rb*rwb+space_back
    for j=1:N
        p(i,j) = as*(0.17*xdim - x(i));
    end
end
for i=(space_front+1.5*x_exit+2*rb*rwb+space_back+1):(2*space_front+2*x_exit+2*rb*rwb+space_back+re1*rwe)
    for j=1:N
        p(i,j) = as*(x(i)-0.2*xdim);
    end
end
for i=(2*space_front+2*x_exit+2*rb*rwb+space_back+re1*rwe+1):(2*space_front+2.5*x_exit+2*rb*rwb+2*space_back+2*re1*rwe)
    for j=1:N
      p(i,j) = as*(0.6*xdim - x(i));   
    end
end
for i=(2*space_front+2.5*x_exit+2*rb*rwb+2*space_back+2*re1*rwe+1):3*space_front+3*x_exit+2*rb*rwb+2*space_back+2*re1*rwe+re2*rwe
    for j=1:N
        p(i,j) = as*(x(i)-0.67*xdim);
    end
end
for i=(3*space_front+3*x_exit+2*rb*rwb+2*space_back+2*re1*rwe+re2*rwe+1):M
    for j=1:N
      p(i,j) = as*(0.93*xdim - x(i));   
    end
end




% row function (modified aisle width) 
%business seating area
for i=(space_front+x_exit+lrb+1):(x_exit+space_front+2*rb*rwb)
    for j=1:NCb
      p(i,j) = p(i,j) + rs*(chairdim_b - ycb(j));
      p(i,N-j+1) = p(i,N-j+1) + rs*(chairdim_b - ycb(j));
    end
    for j=1:((NCb-extension)*3/4+extension) 
      p(i,(N-1)/2.0-j+1) = p(i,(N-1)/2.0-j+1) + rs*(0.87*chairdim_b - ycb(j));
      p(i,(N-1)/2.0+j) = p(i,(N-1)/2.0+j) + rs*(0.87*chairdim_b - ycb(j));
    end
end
%economy cabin1
for i=(2*space_front+space_back+2*x_exit+2*rb*rwb+lre+1):(2*space_front+space_back+2*x_exit+2*rb*rwb+2*re1*rwe)
    for j=1:NCe
      p(i,j) = p(i,j) + rs*(chairdim_e - yce(j));
      p(i,N-j+1) = p(i,N-j+1) + rs*(chairdim_e - yce(j));
    end
    for j=1:((NCe-extension)/2.0 + extension)
      p(i,(N-1)/2.0-j+1) = p(i,(N-1)/2.0-j+1) + rs*(0.6*chairdim_e - yce(j));
      p(i,(N-1)/2.0+j) = p(i,(N-1)/2.0+j) + rs*(0.6*chairdim_e - yce(j));
    end
end
%economy cabin2
for i=(3*space_front+2*space_back+3*x_exit+2*rb*rwb+2*re1*rwe+lre+1):(3*space_front+2*space_back+3*x_exit+2*rb*rwb+2*re1*rwe+2*re2*rwe)
    for j=1:NCe
      p(i,j) = p(i,j) + rs*(chairdim_e - yce(j));
      p(i,N-j+1) = p(i,N-j+1) + rs*(chairdim_e - yce(j));
    end
    for j=1:((NCe-extension)/2.0 + extension)
      p(i,(N-1)/2.0-j+1) = p(i,(N-1)/2.0-j+1) + rs*(0.6*chairdim_e - yce(j));
      p(i,(N-1)/2.0+j) = p(i,(N-1)/2.0+j) + rs*(0.6*chairdim_e - yce(j));
    end
end

%front exit path
for i=1:space_front+x_exit+lrb
    for j=1:(N/2.0)
      p(i,j) = p(i,j) + 1.5*rs*(y(j)-ydim/2.0);
      p(i,N-j+1) = p(i,N-j+1) + 1.5*rs*(y(j)-ydim/2.0);
    end
end
%path between business and 1st economy cabin
for i=(x_exit+space_front+2*rb*rwb+1):(2*space_front+space_back+2*x_exit+2*rb*rwb+lre)
     for j=1:(N/2.0)
      p(i,j) = p(i,j) + 1.5*rs*(y(j)-ydim/2.0);
      p(i,N-j+1) = p(i,N-j+1) + 1.5*rs*(y(j)-ydim/2.0);
    end
end

%path between 1st and 2nd economy cabins
for i=(2*space_front+space_back+2*x_exit+2*rb*rwb+2*re1*rwe+1):(3*space_front+2*space_back+3*x_exit+2*rb*rwb+2*re1*rwe+lre)
    for j=1:(N/2.0)
      p(i,j) = p(i,j) + 1.5*rs*(y(j)-ydim/2.0);
      p(i,N-j+1) = p(i,N-j+1) + 1.5*rs*(y(j)-ydim/2.0);
    end
end

%back exit path
for i=M-x_exit-space_back+1:M
    for j=1:(N/2.0)
      p(i,j) = p(i,j) + 1.5*rs*(y(j)-ydim/2.0);
      p(i,N-j+1) = p(i,N-j+1) + 1.5*rs*(y(j)-ydim/2.0);
    end
end


% seat function
%business class
for i=1:2*rb
    for j=1:ncb
        posx = space_front+x_exit+(i-1)*rwb*numx+1;
        posy = (j-1)*cwb*numy+1;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_b, chair_width_b, cs);
    end
    for j =(ncb+1):3*ncb-1
        posx = space_front+x_exit+(i-1)*rwb*numx+1;
        posy = (j-1)*cwb*numy+ab;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_b, chair_width_b, cs);
    end
    for j =3*ncb:cb-1
        posx = space_front+x_exit+(i-1)*rwb*numx+1;
        posy = (j-1)*cwb*numy+2*ab;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_b, chair_width_b, cs);
    end
end

%1st economy
for i=1:2*re1
    for j=1:nce
        posx = 2*space_front+2*x_exit+space_back+2*rb*rwb+(i-1)*rwe*numx+1;
        posy = (j-1)*cwe*numy+1;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, cs);
    end
    for j =(nce+1):2*nce
        posx = 2*space_front+2*x_exit+space_back+2*rb*rwb+(i-1)*rwe*numx+1;
        posy = (j-1)*cwe*numy+ae;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, cs);
    end
    for j =2*nce+1:3*nce
        posx = 2*space_front+2*x_exit+space_back+2*rb*rwb+(i-1)*rwe*numx+1;
        posy = (j-1)*cwe*numy+2*ae;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, cs);
    end
end
%remove the seat from first row
posx = 2*space_front+2*x_exit+space_back+2*rb*rwb+1;
posy = 2*cwe*numy+1;
p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);
%remove 3 seats from the last row of economy1
posx = 2*space_front+2*x_exit+space_back+2*rb*rwb+(2*re1-1)*rwe*numx+1;
for j = 2*nce+1:3*nce
    posy = (j-1)*cwe*numy+2*ae;
    p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);
end
    


%2nd economy
for i=1:2*re2
    for j=1:nce
        posx = 3*space_front+3*x_exit+2*space_back+2*rb*rwb+2*re1*rwe+(i-1)*rwe*numx+1;
        posy = (j-1)*cwe*numy+1;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, cs);
    end
    for j =(nce+1):2*nce
        posx = 3*space_front+3*x_exit+2*space_back+2*rb*rwb+2*re1*rwe+(i-1)*rwe*numx+1;
        posy = (j-1)*cwe*numy+ae;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, cs);
    end
    for j =2*nce+1:3*nce
        posx = 3*space_front+3*x_exit+2*space_back+2*rb*rwb+2*re1*rwe+(i-1)*rwe*numx+1;
        posy = (j-1)*cwe*numy+2*ae;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, cs);
    end
end
%remove the seats in the front
posx = 3*space_front+3*x_exit+2*space_back+2*rb*rwb+2*re1*rwe+1;
for j=1:nce
    posy = (j-1)*cwe*numy+1;
    p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);
end
for j =2*nce+1:3*nce
    posy = (j-1)*cwe*numy+2*ae;
    p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);
end
%remove the seats in the back
posx = 3*space_front+3*x_exit+2*space_back+2*rb*rwb+2*re1*rwe+(2*re2-2)*rwe*numx+1;
posy = (nce-1)*cwe*numy+1;
p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);

posy = (2*nce)*cwe*numy+2*ae;
p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);

posx = 3*space_front+3*x_exit+2*space_back+2*rb*rwb+2*re1*rwe+(2*re2-1)*rwe*numx+1;
posy = (nce-1)*cwe*numy+1;
p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);

for j = 2*nce+1:3*nce
    posy = (j-1)*cwe*numy+2*ae;
    p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom_e, chair_width_e, -cs);
end

%The very top part-missing row
for i=x_exit:x_exit+space_front+2*rb*rwb+space_back
    p(i,ydim) = cs+p(i,ydim);
end
for i=2*x_exit+space_front+2*rb*rwb+space_back:2*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back
    p(i,ydim) = cs+p(i,ydim);
end
for i=3*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back:M-x_exit
    p(i,ydim) = cs+p(i,ydim);
end
 


%   figure;
%   set(gcf, 'Position', [20 20 3000 900])
%   mesh(y,x,p);
%   axis vis3d
%   pause(2)
%   for i=1:180
%       %clf;
% %       hold on;
%       camorbit(1,0,'data',[0 0 1])
%       drawnow
% %       hold off;
%       set(gcf, 'Position', [20 20 3000 900])
%       %set(gcf,'Position',get(0,'Screensize'));
%       F(i)=getframe;
%   end
  
end