function [ gbest_x, gbest_y ] = gbest( gx,gy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
r = sqrt(gx^2 + gy^2);
gbest_x = gx/(10+r);
gbest_y = gy/(10+r);

end

