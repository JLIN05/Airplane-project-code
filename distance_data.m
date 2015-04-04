function [ r_max,r_ave ] = distance_data( r,r_max,r_count,r_sum );
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%for finding out the range of distance from particle to the nearest exit
      if r>r_max;
          r_max = r;
      end
      
      r_ave=r_sum/r_count;


end

