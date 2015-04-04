function [ grad_count,grad_x_max,grad_x_min,grad_y_max,grad_y_min,grad_ave_x,grad_ave_y,bin_x,bin_y ] = grad_data( pu_x,pu_y,ix,iy,grad_count, grad_x_sum,grad_y_sum,grad_x_max,grad_x_min,grad_y_max,grad_y_min,partition_x, partition_y,bin_x,bin_y )
%for finding out the unscaled grad distribution, and ave after scaled
%   Detailed explanation goes here

%for finding out the range of gradient
            
            if pu_x(ix,iy)>=grad_x_max
                grad_x_max = pu_x(ix,iy);
            end
            if pu_x(ix,iy)<grad_x_min
                grad_x_min = pu_x(ix,iy);
            end
            if pu_y(ix,iy)>=grad_y_max
                grad_y_max = pu_y(ix,iy);
            end
            if pu_y(ix,iy)<grad_y_min
                grad_y_min = pu_y(ix,iy);
            end
            
          %find the distribution of grad
          if pu_x(ix,iy)<partition_x(1)
              bin_x(1,1)=bin_x(1,1)+1;
          end
          if pu_x(ix,iy)>partition_x(length(partition_x))
              bin_x(length(partition_x)+1,1)=bin_x(length(partition_x)+1,1)+1;
          end
          for j = 2:length(partition_x)
              if partition_x(j-1)<=pu_x(ix,iy)&&pu_x(ix,iy)<=partition_x(j)
                bin_x(j,1)=bin_x(j,1)+1;
              end
          end
          
          if pu_y(ix,iy)<partition_y(1)
              bin_y(1,1)=bin_y(1,1)+1;
          end
          if pu_y(ix,iy)>partition_y(length(partition_y))
              bin_y(length(partition_y)+1,1)=bin_y(length(partition_y)+1,1)+1;
          end
          for k = 2:length(partition_y)
              if partition_y(k-1)<=pu_y(ix,iy)&&pu_y(ix,iy)<=partition_y(k)
                bin_y(k,1)=bin_y(k,1)+1;
              end
          end
        
          %for finding out the average of gradient abs
            
            grad_x_sum = grad_x_sum + abs(p_x(ix,iy));
            grad_y_sum = grad_y_sum + abs(p_y(ix,iy));
            
            %finding the average of the grad
           grad_ave_x = grad_x_sum/grad_count;
           grad_ave_y = grad_y_sum/grad_count;

end

