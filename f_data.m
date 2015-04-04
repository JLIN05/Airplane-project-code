function [ f_count,f_max,f_ave,bin_f ] = f_data( f, f_count, f_sum,f_max,partition_f, bin_f )
%for finding out the unscaled pushing forch distribution, and ave after scaled
%   Detailed explanation goes here

%for finding out the range of gradient
            
            if f>=f_max
            f_max = f;
            end
            
          %find the distribution of grad
          if f<partition_f(1)
              bin_f(1,1)=bin_f(1,1)+1;
          end
          if f>partition_f(length(partition_f))
              bin_f(length(partition_f)+1,1)=bin_f(length(partition_f)+1,1)+1;
          end
          for j = 2:length(partition_f)
              if partition_f(j-1)<=f&&f<=partition_f(j)
                bin_f(j,1)=bin_f(j,1)+1;
              end
          end
          
        
          %for finding out the average of gradient abs
            
            f_sum = f_sum + abs(f);
            
            %finding the average of the grad
           f_ave = f_sum/f_count;
           

end

