function [ p,numx,numy,xdim,ydim ] = chair_only(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back)
% stuff to build a plane

% gridspacing
numx = 1;
numy = 1;
dx = 1.0/numx;
dy = 1.0/numy;

%plane dimensions
xdim = r*rw+space_front+space_back; %rows of airplane
ydim = c*cw+a; 
chairdim = nc*cw;

%chair dimensions
seat_depth = sd*numx;
chair_depth = cd*numx;
legroom = lr*numx;
chair_width = cw*numy;

% Set the grid spacing
x = 0:dx:xdim;
y = 0:dy:ydim;
% grid spacing for chair part of plane
yc = 0:dx:chairdim;

% matrix size
M = length(x);
N = length(y);
% chair part of matrix
NC = length(yc);

p = ones(M, N);

% i = row numbers #, x direction
% j = seat numbers #, y direction


% seat function
for i=1:r
    for j=1:nc
        posx = space_front+(i-1)*rw*numx+1;
        posy = (j-1)*cw*numy+1;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width, cs);
    end
    for j =(nc+1):c
        posx = space_front+(i-1)*rw*numx+1;
        posy = (j-1)*cw*numy+a;
        p = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width, cs);
    end
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
for i=1:r
    posx = space_front+(i-1)*rw*numx+1;
    startx = posx + legroom;
    seatx = startx + seat_depth;
    chairbackx = seatx + chair_depth;
    j = ydim+1;
    for k=startx:seatx
       p(k,j) = p(k,j) + 0.25*cs;
    end
    for k=seatx+1:chairbackx
       p(k,j) = p(k,j) + cs;
    end
end


%wall function
for i = 1:wall
    for j = 1:N
        p(i,j) = p(i,j) + ws;
    end
end


%  figure;
%  surf(y,x,p);
%  pause(10.0)



end

