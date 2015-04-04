function [ p,numx,numy,xdim,ydim,x,y ] = block_Boeing737_plane3(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back,space_above,space_below,x_exit,y_exit)
% stuff to build a plane
% 
% gridspacing
numx = 1;
numy = 1;
dx = 1.0/numx;
dy = 1.0/numy;

% % Number of rows
% r = 10;
% % Number of columns
% c = 7; % 3 chairs, 1 for aisle, 3 chairs
% % Aisle width, one chair width
% a = 18;
% 
% % chair width
% cw = 18;
% % row depth
% rw = 32;
% %seat depth
% sd = 18;
% %chair depth
% cd = 4;
% %legroom
% lr = 10;
% 
% 
% % number of chairs per row
% nc = 3;


extension = 9; %for agents to exit faster

%plane dimensions
xdim = 4*r*rw+2*space_front+2*space_back+4*x_exit; %rows of airplane
ydim = c*cw+a+space_above+space_below; 
chairdim = nc*cw+extension;

%chair dimensions
seat_depth = sd*numx;
chair_depth = cd*numx;
legroom = lr*numx;
chair_width = cw*numy;

% Set the grid spacing
x = 0:dx:xdim;
y = 0:dy:ydim;
% grid spacing for chair part of plane
yc = 0:dy:chairdim;


% matrix size
M = length(x);
N = length(y);
% chair part of matrix
NC = length(yc);


p = zeros(M, N);

% i = row numbers #, x direction
% j = seat numbers #, y direction

% aisle function
for i=1:space_front+x_exit+r*rw
    for j=1:N
      p(i,j) = as*(x(i));%/xdim);   
    end
end
for i=(space_front+x_exit+r*rw+1):(M-1)/2.0-x_exit
    for j=1:N
        p(i,j) = as*(0.5*xdim - x(i));
    end
end
for i=((M-1)/2.0-x_exit+1):((M-1)/2.0+x_exit)
    for j=1:N
        p(i,j) = as*(0.5*xdim - x(space_front+x_exit+2*r*rw+space_back));
    end
end
for i=((M-1)/2.0+x_exit+1):(2*space_front+3*x_exit+3*r*rw+space_back)
    for j=1:N
      p(i,j) = as*(x(i)-0.5*xdim);   
    end
end
for i=(2*space_front+3*x_exit+3*r*rw+space_back+1):M
    for j=1:N
        p(i,j) = as*(x(i)-0.5*xdim);
    end
end


% row function (modified aisle width) extention = 6
%front path to exit 1 (tested)
for i=1:space_front+x_exit+lr
    for j=1:NC
      p(i,j) = p(i,j) + 3*rs*(yc(j)-chairdim);
      p(i,N-j+1) = p(i,N-j+1) + 3*rs*(chairdim-yc(j));
    end
end
%first group of passengers
for i=(space_front+x_exit+lr+1):((M-1)/2.0-space_back-x_exit)
    for j=1:NC
      p(i,j) = p(i,j) + rs*(chairdim - yc(j));
      p(i,N-j+1) = p(i,N-j+1) + rs*(chairdim - yc(j));
    end
end
%middle path to exit 2
for i=((M-1)/2.0-space_back-x_exit+1):(M-1)/2.0
    for j=1:(N/2.0)
      p(i,j) = p(i,j) + 3*rs*(y(j)-ydim/2.0);
      p(i,N-j+1) = p(i,N-j+1) + 3*rs*(y(j)-ydim/2.0);
    end
end
%middle path to exit 3
for i=((M-1)/2.0+1):(M-1)/2.0+space_front+x_exit+lr
    for j=1:(N/2.0)
      p(i,j) = p(i,j) + 3*rs*(y(j)-ydim/2.0);
      p(i,N-j+1) = p(i,N-j+1) + 3*rs*(y(j)-ydim/2.0);
    end
end
%second group of passengers
for i=(M-1)/2.0+space_front+x_exit+lr+1:M-x_exit-space_back
    for j=1:NC
      p(i,j) = p(i,j) + rs*(chairdim - yc(j));
      p(i,N-j+1) = p(i,N-j+1) + rs*(chairdim - yc(j));
    end
end

%path to exit 4
for i=M-x_exit-space_back+1:M
    for j=1:NC
      p(i,j) = p(i,j) + 3*rs*(chairdim-yc(j));
      p(i,N-j+1) = p(i,N-j+1) + 3*rs*(chairdim-yc(j));
    end
end


% seat function
for i=1:2*r
    for j=1:nc
        posx = space_front+x_exit+(i-1)*rw*numx+1;
        posy = (j-1)*cw*numy+1;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width, cs);
    end
    for j =(nc+1):c
        posx = space_front+x_exit+(i-1)*rw*numx+1;
        posy = (j-1)*cw*numy+a;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width, cs);
    end
end
for i=2*r+1:4*r
    for j=1:nc
        posx = 2*space_front+3*x_exit+space_back+(i-1)*rw*numx+1;
        posy = (j-1)*cw*numy+1;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width, cs);
    end
    for j =(nc+1):c
        posx = 2*space_front+3*x_exit+space_back+(i-1)*rw*numx+1;
        posy = (j-1)*cw*numy+a;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width, cs);
    end
end
%the middle four seats
 p = add_chair(p, (M-1)/2.0-x_exit, cw*numy+1, seat_depth, chair_depth, legroom, chair_width, cs);
 p = add_chair(p, (M-1)/2.0-x_exit, 2*cw*numy+1, seat_depth, chair_depth, legroom, chair_width, cs);
 p = add_chair(p, (M-1)/2.0-x_exit, 3*cw*numy+a+1, seat_depth, chair_depth, legroom, chair_width, cs);
 p = add_chair(p, (M-1)/2.0-x_exit, 4*cw*numy+a+1, seat_depth, chair_depth, legroom, chair_width, cs);

%The very top part-missing row
for i=x_exit:(M-1)/2.0-x_exit*2
    p(i,ydim) = cs+p(i,ydim);
end
for i=(M-1)/2.0+x_exit*2:M-x_exit
    p(i,ydim) = cs+p(i,ydim);
end
 
% Last chair before aisle
% for i=1:r
%     posx = space_front+(i-1)*rw*numx+1;
%     startx = posx + legroom;
%     seatx = startx + seat_depth;
%     chairbackx = seatx + chair_depth;
%     j = cw*nc+1;
%     for k=startx:seatx
%        p(k,j) = p(k,j) + 0.25*cs;
%     end
%     for k=seatx+1:chairbackx
%        p(k,j) = p(k,j) + cs;
%     end
% end


% Last row of matrix
% for i=1:4*r
%     posx = space_front+(i-1)*rw*numx+1;
%     startx = posx + legroom;
%     seatx = startx + seat_depth;
%     chairbackx = seatx + chair_depth;
%     j = ydim+1;
%     for k=startx:seatx
%        p(k,j) = p(k,j) + 0.25*cs;
%     end
%     for k=seatx+1:chairbackx
%        p(k,j) = p(k,j) + cs;
%     end
% end


% %source of danger
% for i = xdim*0.5+0.5*x_exit:xdim*0.5+1.5*x_exit
%     for j = ydim-space_above-y_exit:ydim-space_above
%         p(i,j) = p(i,j) + 80;
%     end
% end

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

