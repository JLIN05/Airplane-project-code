function [ v_count,v_max,v_ave,bin_v ] = v_data( v, v_count, v_sum,v_max,partition_v, bin_v )
%for finding out the velocities distribution, and ave after scaled
%   Detailed explanation goes here

%for finding out the range of velocities
            if v>=v_max
            v_max = v;
            end
            
          %find the distribution of velocities
          if v<partition_v(1)
              bin_v(1,1)=bin_v(1,1)+1;
          end
          if v>partition_v(length(partition_v))
              bin_v(length(partition_v)+1,1)=bin_v(length(partition_v)+1,1)+1;
          end
          for j = 2:length(partition_v)
              if partition_v(j-1)<=v&&v<=partition_v(j)
                bin_v(j,1)=bin_v(j,1)+1;
              end
          end
          
        
          %for finding out the average of velocity abs
            
            v_sum = v_sum + abs(v);
            
            %finding the average of the grad
           v_ave = v_sum/v_count;
           

end

