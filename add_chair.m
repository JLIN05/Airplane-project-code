function [ p ] = add_chair(p, posx, posy, seat_depth, chair_depth, legroom, chair_width,cs)
% Add a single chair

startx = posx + legroom;
seatx = startx + seat_depth;
chairbackx = seatx + chair_depth;

starty = posy;
endy = starty + chair_width-1;

for j=starty:endy
    for i=startx:seatx
       p(i,j) = p(i,j) + 0.30*cs;
    end
    for i=seatx+1:chairbackx
        p(i,j) = p(i,j) + cs;
    end
end


