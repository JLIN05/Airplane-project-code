% Particle Swarm Optimization Asiana Airplane Simulation 
% Simulates the movements of a swarm to minimize any functions with a 
% moving global minima
%  
% The swarm matrix is  
%
% swarm(index, [1-location, 2-velocity, 3-emotion, 4-best 
% value], [1-x or the value component, 2-y component])

% Initialize random number generator
random_number_generator = rng('shuffle');

% Initialization
% Parameters
iterations = 8000;
w = 0.85; %inertia factor
c1 = 1.5;%2.5; %self-confidence factor
c2 = 0.5;%1.0; %swarm confidence
swarm_size = 172;

f = zeros(swarm_size,1);


dt = 0.05;%0.25; %time-step
% Number of rows
r = 7;
% Number of columns
c = 6; % 3 chairs, (1 for aisle), 3 chairs
% Aisle width, one chair width
a = 18;


% chair width
cw = 18;
% row depth (lr+sd+cd)
rw = 50;%32;
%seat depth
sd = 18;
%chair depth
cd = 11; %4
%legroom
lr = 21;
%thickness of the wall
wall = 0;%15;%15;
%space in the front of the plane
space_front = 25;%2000;%40;
%space in the back of the plane
space_back = 46;%2000;%40;
%space for the particles to run off the screen
space_out = 0;%1960;%20;
%space above
space_above = 0;%500;
%space below
space_below = 0;%500;
%width of the exits
x_exit = 40;%20;
%length of the exits
y_exit = 10;
%person radius
rad = 9;
% personal comfort tolerance
tol = 2*rad + 6.0;

% aisle stength
as = 1.0;%1.0;60; %20;
% chair strength
cs = 380;%130;%60;%30;%20;
% row stength 
rs = 1.0;%10;
%height of the wall
ws = 30;

% repulsion strength
c3 = 1;%80.0;%40;35;%65

% number of chairs per row
nc = 3;

%set the potential function of the plane
[ p,numx,numy,xdim,ydim ]=Boeing737_plane(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back,space_above,space_below,x_exit,y_exit);
%adding chair
[ pc]=chair_Boeing737(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back,space_above,space_below,x_exit);
%gradiant for particle local search
[p_x,p_y] = plane_grad( p,rad,numx,numy );
%[pu_x,pu_y] = plane_grad_unscaled( p,rad,numx,numy );

%for finding out range of gradients
grad_x_max = 0;
grad_y_max = 0;
grad_x_min = grad_x_max;
grad_y_min = grad_y_max;

%for finding out average of gradients
grad_count = 0;
grad_x_sum = 0;
grad_y_sum = 0;

%finding the distribution of x_grad and y_grad
partition_x = -4:0.5:4;
partition_y = -2:0.5:2;
bin_x = zeros(length(partition_x)+1,1);
bin_y = zeros(length(partition_y)+1,1);

%for finding out the pushing force
f_max = -10;
f_sum = 0;
f_count = 0;
partition_f = 0: 0.1:0.7;
bin_f = zeros(length(partition_f)+1,1);

%for finding out average of r
r_max = 0;
r_sum = 0;
r_count = 0;

%for finding out the particles velocies
v_max = 0;
v_sum = 0;
v_count = 0;
partition_v = 0:7;
bin_v = zeros(length(partition_v)+1,1);

%preset the matrix
swarm = zeros(swarm_size, 4, 2);
s = zeros(swarm_size, 4, 2);

%for plotting the seats
x1_s = (lr+sd+space_front+x_exit): rw : rw*2*r+space_front+x_exit;
x1_e = rw+space_front+x_exit: rw: rw*2*r+space_front+x_exit;
x2_s = (lr+sd+rw*2*r+2*space_front+3*x_exit+space_back): rw : rw*4*r+2*space_front+3*x_exit+space_back;
x2_e = rw+rw*2*r+2*space_front+3*x_exit+space_back: rw: rw*4*r+2*space_front+3*x_exit+space_back;
y_s = space_below: (cw*nc+a) : (nc*cw+a+space_below);
y_e = cw*nc+space_below: (cw*nc+a) : (nc*cw+a+space_below);

% initial swarm position and velocity

    index = 1;
    for i = 1 : (2*r)
        for j = 1 : 2*nc
            if ceil(j/nc)==1
            swarm(index, 1, 1) = space_front+x_exit+rw*i-cd-sd/2;
            swarm(index, 1, 2) = -cw/2+cw*j+space_below;
            elseif ceil(j/nc)==2
            swarm(index, 1, 1) = space_front+x_exit+rw*i-cd-sd/2;
            swarm(index, 1, 2) = -cw/2+a+cw*j+space_below;
            end
            index = index + 1;
            if(index>swarm_size) 
                break; end
        end  
        if(index>swarm_size) 
            break; end
    end
if swarm_size>(172-4)/2
    index = (172-4)/2+1;
    for i = (2*r+1) : (4*r)
        for j = 1 : 2*nc
            if ceil(j/nc)==1
            swarm(index, 1, 1) = 2*space_front+3*x_exit+space_back+rw*i-cd-sd/2;
            swarm(index, 1, 2) = -cw/2+cw*j+space_below;
            elseif ceil(j/nc)==2
            swarm(index, 1, 1) = 2*space_front+3*x_exit+space_back+rw*i-cd-sd/2;
            swarm(index, 1, 2) = -cw/2+a+cw*j+space_below;
            end
            index = index + 1;
            if(index>swarm_size) 
                break; end
        end  
        if(index>swarm_size) 
            break; end
    end
end
%putting 4 ppl in the middle

for i = 0: 1
swarm(swarm_size-i, 1, 1) = space_front+space_back+x_exit+rw*(2*r+1)-cd-sd/2;
swarm(swarm_size-i, 1, 2) = -cw/2+(nc+i-1)*cw+space_below;
end

for i = 2: 3
swarm(swarm_size-i, 1, 1) = space_front+space_back+x_exit+rw*(2*r+1)-cd-sd/2;
swarm(swarm_size-i, 1, 2) = -cw/2+(nc+i-1)*cw+a+0.5+space_below;
end



swarm(:, 4, 1) = 1000;          % best value so far
swarm(:, 2, :) = 0;             % initial velocity

% Atttraction length
La = 48;%36
% Repulsion length
Lr = 24;%18
% Attraction strength
Ca = 5; %30
% Repultion strength
Cr = 15;%50



% %outputting video
%  writerObj = VideoWriter('737_no emo.avi');
%  writerObj = VideoWriter('C:\Users\Junyuan Lin\Dropbox\Pepp\plane project_new\737_no emo', 'MPEG-4');
% % %witerObj = VideoWriter('plane.avi','Uncompressed AVI');
%  writerObj.FrameRate = 4;
%  writerObj.Quality= 100;
%  open(writerObj)

   
% Here make a matrix of all global best positions (8 in this case)
% call that matrix exits and then have each row be an x,y position of
% an exit
exits = zeros(8,2);
exits(1,1) = space_out+wall+x_exit;
exits(2,1) = space_out+wall+x_exit;
exits(3,1) = xdim-space_out-wall-x_exit;
exits(4,1) = xdim-space_out-wall-x_exit;
exits(5,1) = xdim/2.0-x_exit;
exits(6,1) = xdim/2.0-x_exit;
exits(7,1) = xdim/2.0+x_exit;
exits(8,1) = xdim/2.0+x_exit;
exits(1,2) = space_below;%+y_exit/2;
exits(2,2) = ydim-space_above;%-y_exit/2;
exits(3,2) = space_below;%+y_exit/2;
exits(4,2) = ydim-space_above;%-y_exit/2;
exits(5,2) = space_below;%+y_exit/2;
exits(6,2) = ydim-space_above;%-y_exit/2;
exits(7,2) = space_below;%+y_exit/2;
exits(8,2) = ydim-space_above;%-y_exit/2;


dis = zeros(8,1);
xdis = zeros(8,1);
ydis = zeros(8,1);




figure
hold on

pfig=plot(1);
qfig=plot(1);

%set seats and chairs
for j = 1:length(x1_s)
    for k = 1:length(y_s)
        for m = 1:nc
            rectangle('Position',[x1_s(j),y_s(k),cd,cw*m]);
            rectangle('Position',[x1_s(j)-sd,y_s(k)+(m-1)*cw,sd,cw],'LineWidth',2,'LineStyle','--');
        end
    end
end
for j = 1:length(x2_s)
    for k = 1:length(y_s)
        for m = 1:nc
            rectangle('Position',[x2_s(j),y_s(k),cd,cw*m]);
            rectangle('Position',[x2_s(j)-sd,y_s(k)+(m-1)*cw,sd,cw],'LineWidth',2,'LineStyle','--');
        end
    end
end

%the middle four seats and chairs
for m = 1:nc-1
    rectangle('Position',[xdim/2.0,cw,cd,m*cw]);
    rectangle('Position',[xdim/2.0-sd,cw*m,sd,cw],'LineWidth',2,'LineStyle','--');
    
    rectangle('Position',[xdim/2.0,nc*cw+a,cd,m*cw]);
    rectangle('Position',[xdim/2.0-sd,nc*cw+a+(m-1)*cw,sd,cw],'LineWidth',2,'LineStyle','--');
end

%plotting exits
%rectangle('Position',[0,0,wall,ydim],'FaceColor','r')
rectangle('Position',[space_out+wall,space_below,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[space_out+wall,ydim-space_above-y_exit,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim-(space_out+x_exit+wall),space_below,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim-(space_out+x_exit+wall),ydim-space_above-y_exit,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim/2.0-2*x_exit,space_below,2*x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim/2.0-2*x_exit,ydim-space_above-y_exit,2*x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim/2.0,space_below,2*x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim/2.0,ydim-space_above-y_exit,2*x_exit,y_exit],'FaceColor','g')
axis([space_out+wall xdim-space_out space_below ydim-space_above]);
%grid off
set(gcf, 'Position', [100 100 2.5*(xdim-space_out) 5*(ydim-space_above-space_below)]);
frame = getframe;
fig = gcf;







% Iterations
iter = 1;
for iter = 1: iterations
    fx = zeros(swarm_size,1);
    fy = zeros(swarm_size,1); 
    
    for i = 1 : swarm_size
        
        s(i, 1, 1) = swarm(i, 1, 1) + 6.5*swarm(i, 2, 1)*dt;     %update x position
        s(i, 1, 2) = swarm(i, 1, 2) + 6.5*swarm(i, 2, 2)*dt;     %update y position
%         
%       prevent agents to go off screen
        if s(i, 1, 2) >=ydim-rad || s(i, 1, 2) <=rad
            swarm(i, 2, 2) = -0.2*swarm(i, 2, 2);
        end
        if (s(i, 1, 1))<=rad||s(i, 1, 1)>=xdim-rad
            swarm(i, 2, 1) = -0.2*swarm(i, 2, 1);
        end
     
     end
    

    for i=1:swarm_size
        
        swarm(i, 1, 1) = swarm(i, 1, 1) + 6.5*swarm(i, 2, 1)*dt;     %update x position
        swarm(i, 1, 2) = swarm(i, 1, 2) + 6.5*swarm(i, 2, 2)*dt;     %update y position
            
        x = swarm(i, 1, 1);
        y = swarm(i, 1, 2);
        
        for j = (i+1) : swarm_size
                dx = swarm(j,1,1)-swarm(i,1,1);
                dy = swarm(j,1,2)-swarm(i,1,2);
                
                r = sqrt(dx*dx + dy*dy);
                
                if( r < tol)
                
                    x2 = round((x+dx/2)*numx);
                    y2 = round((y+dy/2)*numy);
                    chair_repel = pc(x2,y2);
                
                    u = c3*(Ca/La * exp(-r / La) - Cr/Lr * exp(-r / Lr));
                    %u = -c3/r^2/chair_repel;
                    fx(i) = fx(i) + u*dx/r;
                    fy(i) = fy(i) + u*dy/r;
                    fx(j) = fx(j) - u*dx/r;
                    fy(j) = fy(j) - u*dy/r; 
                end
         end
 
%         f(i) = sqrt(fx(i)*fx(i)+fy(i)*fy(i));
%         f_count = f_count+1;
%         [ f_count,f_max,f_ave,bin_f ] = f_data( f(i), f_count, f_sum,f_max,partition_f, bin_f );
        
        
        ix=round(x*numx);
        iy=round(y*numy);
        swarm(i,4,1)=p(ix,iy); % fitness evaluation

        
    end
   
    
    deletions = [];
    for i = 1 : swarm_size 
        
          x = swarm(i, 1, 1);
          y = swarm(i, 1, 2);
  
          ix = round(x*numx);
          iy = round(y*numy);
          
          
            
%           [ grad_count,grad_x_max,grad_x_min,grad_y_max,grad_y_min,grad_ave_x,grad_ave_y,bin_x,bin_y ] = grad_data( ix,iy,grad_count, grad_x_sum,grad_y_sum,grad_x_max,grad_x_min,grad_y_max,grad_y_min,partition_x, partition_y,bin_x,bin_y );
%             grad_count = grad_count + 1;
            
            
          
          
          % here compute the distance to each global best.
          % then choose the gbest that is the shortest distance
          % then call that gbestpos(suggesting exit)
          for j=1:8  
            xdis(j) = exits(j,1)-swarm(i,1,1);
            ydis(j) = exits(j,2)-swarm(i,1,2);
            dis(j) = xdis(j)*xdis(j) + ydis(j)*ydis(j);
          end
          [r2,pos]=min(dis);
          r = sqrt(r2);
          
          % if at the exit, add index to a vector called deletions
         
          if r<=2.8*rad
              deletions = i;
          end
          
%           r_count = r_count+1;
%           r_sum = r_sum+r;
          %[ r_max,r_ave ] = distance_data( r,r_max,r_count,r_sum );
          
          
          
          % This is more efficient than gbest function.  
          gbest_x = xdis(pos)/(100.0+r);
          gbest_y = ydis(pos)/(100.0+r);
          
          ux1 = rand();
          ux2 = rand();
          uy1 = rand();
          uy2 = rand();
          
          swarm(i, 2, 1) = fx(i) + w*swarm(i,2,1) - c1*ux1*p_x(ix,iy) + c2*ux2*gbest_x;   
          swarm(i, 2, 2) = fy(i) + w*swarm(i,2,2) - c1*uy1*p_y(ix,iy) + c2*uy2*gbest_y;
          
          v(i) = sqrt(swarm(i,2,1)*swarm(i,2,1)+swarm(i,2,2)*swarm(i,2,2));
          v_count = v_count+1;
          [ v_count,v_max,v_ave,bin_v ] = v_data( v(i), v_count, v_sum,v_max,partition_v, bin_v );
        
    end
        %deleting the agent who exits the plane
        swarm(deletions,:,:)=[];
        swarm_size=length(swarm(:,1,1));
        if swarm_size == 0
            total_iter = iter;
        end
        
    
    % Plotting the swarm
    
    if mod(iter,1)==0
%        clf
%        hold on
%         for i = 1: swarm_size
%             scatter(swarm(i, 1, 1),swarm(i, 1, 2),800,[f(i) 0 1],'filled');
%         end

        delete(pfig);
        delete(qfig);

        pfig = plot(swarm(:, 1, 1), swarm(:, 1, 2), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','MarkerSize',4.0*rad); % drawing swarm movements
        qfig = quiver(swarm(:,1,1),swarm(:,1,2),swarm(:,2,1),swarm(:,2,2),0.5,'k');

         grid off
         set(gcf, 'Position', [100 100 2.5*(xdim-space_out) 5*(ydim-space_above-space_below)]);
         frame = getframe;
         %writeVideo(writerObj,frame);
        
        %pause(4)
    end
    if(swarm_size == 0) 
            break; end
end

hold off




% close(writerObj);
%% taken out the video output part
 