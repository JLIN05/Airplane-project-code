% Particle Swarm Optimization Asiana Airplane Simulation 
% Simulates the movements of a swarm to minimize any functions with a 
% moving global minima
%  
% The swarm matrix is  
%
% swarm(index, [1-location, 2-velocity, 3-emotion, 4-best 
% value], [1-x or the value component, 2-y component])

% Initialization
% Parameters
iterations = 1000;
w = 0.85; %inertia factor
c1 = 2.0; %self-confidence factor
c2 = 1; %swarm confidence
swarm_size = 60;


dt = 0.25; %time-step
% Number of rows
r = 10;
% Number of columns
c = 6; % 3 chairs, (1 for aisle), 3 chairs
% Aisle width, one chair width
a = 18;


% chair width
cw = 18;
% row depth
rw = 32;
%seat depth
sd = 0; %18
%chair depth
cd = 11; %4
%legroom
lr = 21;
%thickness of the wall
wall = 15;
%space in the front of the plane
space_front = 21;
%space in the back of the plane
space_back = 850;%;21
%person radius
rad = 11;

% aisle stength
as = 20;%1.0;
% chair strength
cs = 20;
% row stength 
rs = 10;%1.0;
%height of the wall
ws = 30;

% repulsion strength
repel = 30;%35;

% number of chairs per row
nc = 3;

[ p,numx,numy,xdim,ydim ]=plane(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back);
[ pc]=chair_only(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back);
[p_x,p_y] = plane_grad( p,rad,numx,numy );

%for finding out average of gradients
grad_count = 0;
grad_x_sum = 0;
grad_y_sum = 0;


swarm = zeros(swarm_size, 4, 2);
s = zeros(swarm_size, 4, 2);
x_s = (lr+sd+space_front): rw : rw*r+space_front;
x_e = rw+space_front: rw: rw*r+space_front;
y_s = 0: (cw*nc+a) : (nc*cw+a);
y_e = cw*nc: (cw*nc+a) : (nc*cw+a);

% initial swarm position and velocity
index = 1;
for i = 1 : r
    for j = 1 : 2*nc
        if ceil(j/nc)==1
        swarm(index, 1, 1) = space_front+rw*i-cd-sd/2;
        swarm(index, 1, 2) = -cw/2+cw*j;
        elseif ceil(j/nc)==2
        swarm(index, 1, 1) = space_front+rw*i-cd-sd/2;
        swarm(index, 1, 2) = -cw/2+a+cw*j;
        end
        index = index + 1;
        if(index>swarm_size) 
            break; end
    end  
    if(index>swarm_size) 
        break; end
end


swarm(:, 4, 1) = 1000;          % best value so far
swarm(:, 2, :) = 0;             % initial velocity

% Atttraction length
La = 40;%36
% Repultion length
Lr = 20;%18
% Attraction strength
Ca = 20; %30
% Repultion strength
Cr = 50;%50

% Iterations

%  writerObj = VideoWriter('plane.avi');
%  writerObj = VideoWriter('C:\Users\Junyuan Lin\Dropbox\Pepp\plane project\plane', 'MPEG-4');
% %witerObj = VideoWriter('plane.avi','Uncompressed AVI');
%  writerObj.FrameRate = 4;
%  writerObj.Quality= 100;
%  open(writerObj)

iter = 1;
while min(swarm(:, 1, 1))<(xdim-790)
    fx = zeros(swarm_size,1);
    fy = zeros(swarm_size,1); 
    
    for i = 1 : swarm_size
        
        s(i, 1, 1) = swarm(i, 1, 1) + swarm(i, 2, 1)*dt;     %update x position
        s(i, 1, 2) = swarm(i, 1, 2) + swarm(i, 2, 2)*dt;     %update y position
        
        if  s(i, 1, 1)>=xdim-rad      
            swarm(i, 2, 1) = 0;           
            swarm(i, 2, 2) = 0;
        end
        
        swarm(i, 1, 1) = swarm(i, 1, 1) + swarm(i, 2, 1)*dt;     %update x position
        swarm(i, 1, 2) = swarm(i, 1, 2) + swarm(i, 2, 2)*dt;     %update y position
              
        
        
        x = swarm(i, 1, 1);
        y = swarm(i, 1, 2);
        
        for j = (i+1) : swarm_size
                dx = swarm(j,1,1)-swarm(i,1,1);
                dy = swarm(j,1,2)-swarm(i,1,2);
                
                r = sqrt(dx*dx + dy*dy);
                x2 = round((x+dx/2)*numx);
                y2 = round((y+dy/2)*numy);
                chair_repel = pc(x2,y2);
                
                %u = repel*(Ca/La * exp(-r / La) - Cr/Lr * exp(-r / Lr));
                u = -repel/r^2/chair_repel;
                fx(i) = fx(i) + u*dx/r;
                fy(i) = fy(i) + u*dy/r;
                fx(j) = fx(j) - u*dx/r;
                fy(j) = fy(j) - u*dy/r; 
                
         end
        
      
        ix=round(x*numx);
        iy=round(y*numy);
        swarm(i,4,1)=p(ix,iy); % fitness evaluation
        
        %for finding out the average of gradient abs
            
%             grad_count = grad_count + 1;
%             grad_x_sum = grad_x_sum + abs(p_x(ix,iy));
%             grad_y_sum = grad_y_sum + abs(p_y(ix,iy));
            
     
       
    end

    [temp, gbestpos] = min(swarm(:, 4, 1));        % global best position
    
   
    for i = 1 : swarm_size 
        
          x = swarm(i, 1, 1);
          y = swarm(i, 1, 2);
          
          ix = round(x*numx);
          iy = round(y*numy);
          
          gx(i)=swarm(gbestpos,1,1)-swarm(i,1,1);
          gy(i)=swarm(gbestpos,1,2)-swarm(i,1,2); 
          [gbest_x, gbest_y] = gbest(gx(i),gy(i));
          
          ux1 = rand();
          ux2 = rand();
          uy1 = rand();
          uy2 = rand();
          
          swarm(i, 2, 1) = fx(i) + w*swarm(i,2,1) - c1*ux1*p_x(ix,iy) + c2*ux2*gbest_x;   
          swarm(i, 2, 2) = fy(i) + w*swarm(i,2,2) - c1*uy1*p_y(ix,iy) + c2*uy2*gbest_y;

    end
    

   
    % Plotting the swarm
    if mod(iter,5)==0
        clf
        hold on
        plot(swarm(:, 1, 1), swarm(:, 1, 2), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','MarkerSize',3*rad); % drawing swarm movements
        quiver(swarm(:,1,1),swarm(:,1,2),swarm(:,2,1),swarm(:,2,2),0.5);
        for j = 1:length(x_s)
            for k = 1:length(y_s)
                rectangle('Position',[x_s(j),y_s(k),cd,nc*cw]);
            end
        end
        rectangle('Position',[0,0,wall,ydim],'FaceColor','r')
        rectangle('Position',[xdim-810,nc*cw,20,a],'FaceColor','g')
        axis([0 xdim-790 0 ydim]);
        grid off
        set(gcf, 'Position', [100 100 2.5*(xdim-790) 5*ydim]);
        frame = getframe;
        %writeVideo(writerObj,frame);
        
        %pause(4)
    end
iter=iter+1;
end

%close(writerObj);
%% taken out the video output part
 