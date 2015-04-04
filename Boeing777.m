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
rand(1,100);

% Initialization
% Parameters
iterations = 4000;
w = 0.85; %inertia factor
c1 = 1.5;%2.5; %self-confidence factor
c2 = 0.2; %swarm confidence
swarm_size = 19+158+114;

f = zeros(swarm_size,1);


dt = 0.05; %time-step

% data for business cabin
% Number of rows
rb = 2;
% Number of columns
cb = 8; % 2 chairs, (1 for aisle), 1 chair + 1 half-length chair, 1 half-length chair+ 1 chair, (1 for aisle) 2 chairs
% unit number of chairs per row
ncb = 2;
% Aisle width (to match up with the ydim of economy cabin)
ab = 29;
% chair width
cwb = 20;
%seat depth
sd = 18;
%chair depth
cd = 11; %4
%legroom
lrb = 62;
% row depth (lr+sd+cd)
rwb = lrb+sd+cd;


% data for economy cabin1
% Number of rows
re1 = 9;
% Number of columns
ce = 12; % 3 chairs, (1 for aisle), 3 half-length chairs, 3 half-length chairs, (1 for aisle) 3 chairs
% unit number of chairs per row
nce = 3;
% Aisle width, one chair width
ae = 18;
% chair width
cwe = 18;
%legroom
lre = 33;
% row depth (lr+sd+cd)
rwe = lre+sd+cd;

% data for economy cabin2
% Number of rows
re2 = 7;


%space in the front of the plane
space_front = 25;%2000;%40;
%space in the back of the plane
space_back = space_front+lrb;%2000;%40;
%space for the particles to run off the screen
space_out = 0;%1960;%20;
%space above
space_above = 0;%500;
%space below
space_below = 0;%500;
%width of the exits
x_exit = 80;%20;
%length of the exits
y_exit = 10;
%for deleting agents
deletions = [];
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

% repulsion strength
c3 = 1;%80.0;%40;35;%65


%set the potential function of the plane
[ p,numx,numy,xdim,ydim ]=Boeing777_plane_block_one3(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit);
%adding chair
[ pc]=chair_Boeing777(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit);
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
xb_s = (lrb+sd+space_front+x_exit): rwb : rwb*2*rb+space_front+x_exit;
xe1_s = (lre+sd+rwb*2*rb+2*space_front+2*x_exit+space_back): rwe : rwe*2*re1+rwb*2*rb+2*space_front+2*x_exit+space_back;
xe2_s = (lre+sd+rwe*2*re1+rwb*2*rb+3*space_front+3*x_exit+2*space_back): rwe : xdim-x_exit-space_back;

yb_s_l = [0 cwb];
yb_s_m = [ncb*cwb+ab (ncb+1)*cwb+ab (ncb+2)*cwb+ab];
yb_s_r = [(ncb+3)*cwb+2*ab (ncb+4)*cwb+2*ab];

ye_s_l2 = [0 cwe];
ye_s_l = [0 cwe 2*cwe];
ye_s_m = [nce*cwe+ae nce*cwe+ae+cwe nce*cwe+ae+2*cwe];
ye_s_r = [2*nce*cwe+2*ae 2*nce*cwe+2*ae+cwe 2*nce*cwe+2*ae+2*cwe];
ye_s_r2 = [2*nce*cwe+2*ae 2*nce*cwe+2*ae+cwe];

% initial swarm position and velocity
index = 1;
for i = 1 : 2*rb
    for j = 1 : 2*ncb
        if ceil(j/ncb)==1
        swarm(index, 1, 1) = space_front+x_exit+rwb*i-cd-sd/2;
        swarm(index, 1, 2) = -cwb/2+cwb*j;
        elseif ceil(j/ncb)==2
        swarm(index, 1, 1) = space_front+x_exit+rwb*i-cd-sd/2;
        swarm(index, 1, 2) = -cwb/2+2*ab+3*cwb+cwb*j;
        end
        index = index + 1;
        if(index>swarm_size) 
            break; end
    end  
    if(index>swarm_size)
        break; end
end
if swarm_size>=2*rb*2*ncb+3
        for k = 1:3
            swarm(index,1,1) = space_front+x_exit+rwb-cd-sd/2;
            swarm(index,1,2) = -cwb/2+ab+ncb*cwb+cwb*k;
            index = index + 1;
        end
end
if swarm_size>=2*rb*2*ncb+3+158
    %first row in economy 1
    for i = 1:2
        swarm(index, 1, 1) = 2*space_front+2*x_exit+2*rb*rwb+space_back+rwe-cd-sd/2;
        swarm(index, 1, 2) = -cwe/2+cwe*i;
        index = index + 1;
    end
    for i = 1:3
        swarm(index, 1, 1) = 2*space_front+2*x_exit+2*rb*rwb+space_back+rwe-cd-sd/2;
        swarm(index, 1, 2) = -cwe/2+cwe*nce+ae+cwe*i;
        index = index + 1;
    end
    for i = 1:3
        swarm(index, 1, 1) = 2*space_front+2*x_exit+2*rb*rwb+space_back+rwe-cd-sd/2;
        swarm(index, 1, 2) = -cwe/2+2*cwe*nce+2*ae+cwe*i;
        index = index + 1;
    end
    %row 2-17 in economy 1
    for i = 2 : (2*re1-1)
        for j = 1 : 3*nce
            if ceil(j/nce)==1
            swarm(index, 1, 1) = 2*space_front+2*x_exit+space_back+2*rb*rwb+rwe*i-cd-sd/2;
            swarm(index, 1, 2) = -cwe/2+cwe*j;
            elseif ceil(j/nce)==2
            swarm(index, 1, 1) = 2*space_front+2*x_exit+space_back+2*rb*rwb+rwe*i-cd-sd/2;
            swarm(index, 1, 2) = -cwe/2+ae+cwe*j;
            elseif ceil(j/nce)==3
            swarm(index, 1, 1) = 2*space_front+2*x_exit+space_back+2*rb*rwb+rwe*i-cd-sd/2;
            swarm(index, 1, 2) = -cwe/2+2*ae+cwe*j;
            end
            index = index + 1;
            if(index>swarm_size) 
                break; end
        end  
        if(index>swarm_size) 
            break; end
    end
    
    %row 18
    for i = 1:3
       swarm(index, 1, 1) = 2*space_front+2*x_exit+space_back+2*rb*rwb+rwe*2*re1-cd-sd/2; 
       swarm(index, 1, 2) = -cwe/2+cwe*i;
       index = index + 1;
    end
    for i = 4:6
       swarm(index, 1, 1) = 2*space_front+2*x_exit+space_back+2*rb*rwb+rwe*2*re1-cd-sd/2; 
       swarm(index, 1, 2) = -cwe/2+ae+cwe*i;
       index = index + 1;
    end
    
end


% economy 2
if swarm_size>=2*rb*2*ncb+3+158+114
    %first row in economy 2
    for i = 1:3
        swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe-cd-sd/2;
        swarm(index, 1, 2) = -cwe/2+cwe*nce+ae+cwe*i;
        index = index + 1;
    end
    
    %row 2-12 in economy 2
    for i = 2 : (2*re2-2)
        for j = 1 : 3*nce
            if ceil(j/nce)==1
            swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*i-cd-sd/2;
            swarm(index, 1, 2) = -cwe/2+cwe*j;
            elseif ceil(j/nce)==2
            swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*i-cd-sd/2;
            swarm(index, 1, 2) = -cwe/2+ae+cwe*j;
            elseif ceil(j/nce)==3
            swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*i-cd-sd/2;
            swarm(index, 1, 2) = -cwe/2+2*ae+cwe*j;
            end
            index = index + 1;
            if(index>swarm_size) 
                break; end
        end  
        if(index>swarm_size) 
            break; end
    end
    %last 2 rows
    for i = 1:2
       swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*(2*re2-1)-cd-sd/2; 
       swarm(index, 1, 2) = -cwe/2+cwe*i;
       index = index + 1;
    end
    for i = 1:3
       swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*(2*re2-1)-cd-sd/2; 
       swarm(index, 1, 2) = -cwe/2+ae+cwe*(nce+i);
       index = index + 1;
    end
    for i = 1:2
        swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*(2*re2-1)-cd-sd/2;
        swarm(index, 1, 2) = -cwe/2+2*ae+cwe*(2*nce+i);
        index = index + 1;
    end
    for i = 1:2
       swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*2*re2-cd-sd/2;
       swarm(index, 1, 2) = -cwe/2+cwe*i;
       index = index + 1;
    end
    for i = 1:3
       swarm(index, 1, 1) = 3*space_front+3*x_exit+2*rb*rwb+2*re1*rwe+2*space_back+rwe*2*re2-cd-sd/2;
       swarm(index, 1, 2) = -cwe/2+ae+cwe*(nce+i);
       index = index + 1;
    end
    
end





swarm(:, 4, 1) = 1000;          % best value so far
swarm(:, 2, :) = 0;             % initial velocity

% Atttraction length
La = 48;%36
% Repultion length
Lr = 24;%18
% Attraction strength
Ca = 5; %30
% Repultion strength
Cr = 15;%50



% %outputting video
%  writerObj = VideoWriter('777_no emo_3.30.avi');
%  writerObj = VideoWriter('C:\Users\Junyuan Lin\Dropbox\Pepp\plane project_new\777_no emo_3.30', 'MPEG-4');
% % %witerObj = VideoWriter('plane.avi','Uncompressed AVI');
%  writerObj.FrameRate = 1;
%  writerObj.Quality= 100;
%  open(writerObj)

   
% Here make a matrix of all global best positions (8 in this case)
% call that matrix exits and then have each row be an x,y position of
% an exit
exits = zeros(8,2);
exits(1,1) = x_exit/2;
exits(2,1) = x_exit/2;
exits(3,1) = xdim-x_exit/2;
exits(4,1) = xdim-x_exit/2;
exits(5,1) = x_exit+space_front+2*rb*rwb+space_back+0.5*x_exit;
exits(6,1) = x_exit+space_front+2*rb*rwb+space_back+0.5*x_exit;
exits(7,1) = 2.5*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back;
exits(8,1) = 2.5*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back;
exits(1,2) = y_exit/2;
exits(2,2) = ydim-y_exit/2;
exits(3,2) = y_exit/2;
exits(4,2) = ydim-y_exit/2;
exits(5,2) = y_exit/2;
exits(6,2) = ydim-y_exit/2;
exits(7,2) = y_exit/2;
exits(8,2) = ydim-y_exit/2;


dis = zeros(8,1);
xdis = zeros(8,1);
ydis = zeros(8,1);





figure
hold on

pfig=plot(1);
qfig=plot(1);


%set seats and chairs
%business class
for j = 1:length(xb_s)
    for k = 1:length(yb_s_l)
        rectangle('Position',[xb_s(j),yb_s_l(k),cd,cwb]);
        rectangle('Position',[xb_s(j)-sd,yb_s_l(k),sd,cwb],'LineWidth',2,'LineStyle','--');
    end
end
for j = 1:length(xb_s)
    for k = 1:length(yb_s_m)
        rectangle('Position',[xb_s(j),yb_s_m(k),cd,cwb]);
        rectangle('Position',[xb_s(j)-sd,yb_s_m(k),sd,cwb],'LineWidth',2,'LineStyle','--');
    end
end
for j = 1:length(xb_s)
    for k = 1:length(yb_s_r)
        rectangle('Position',[xb_s(j),yb_s_r(k),cd,cwb]);
        rectangle('Position',[xb_s(j)-sd,yb_s_r(k),sd,cwb],'LineWidth',2,'LineStyle','--');
    end
end

%economy1
for k = 1:length(ye_s_l2)
    rectangle('Position',[xe1_s(1),ye_s_l2(k),cd,cwe]);
    rectangle('Position',[xe1_s(1)-sd,ye_s_l2(k),sd,cwe],'LineWidth',2,'LineStyle','--');
end
for j = 2:length(xe1_s)
    for k = 1:length(ye_s_l)
        rectangle('Position',[xe1_s(j),ye_s_l(k),cd,cwe]);
        rectangle('Position',[xe1_s(j)-sd,ye_s_l(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end

for j = 1:length(xe1_s)
    for k = 1:length(ye_s_m)
        rectangle('Position',[xe1_s(j),ye_s_m(k),cd,cwe]);
        rectangle('Position',[xe1_s(j)-sd,ye_s_m(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end

for j = 1:length(xe1_s)-1
    for k = 1:length(ye_s_r)
        rectangle('Position',[xe1_s(j),ye_s_r(k),cd,cwe]);
        rectangle('Position',[xe1_s(j)-sd,ye_s_r(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end

%economy2
for j = 2:length(xe2_s)-2
    for k = 1:length(ye_s_l)
        rectangle('Position',[xe2_s(j),ye_s_l(k),cd,cwe]);
        rectangle('Position',[xe2_s(j)-sd,ye_s_l(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end
for j = length(xe2_s)-1:length(xe2_s)
    for k = 1:length(ye_s_l2)
        rectangle('Position',[xe2_s(j),ye_s_l2(k),cd,cwe]);
        rectangle('Position',[xe2_s(j)-sd,ye_s_l2(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end

for j = 1:length(xe2_s)
    for k = 1:length(ye_s_m)
        rectangle('Position',[xe2_s(j),ye_s_m(k),cd,cwe]);
        rectangle('Position',[xe2_s(j)-sd,ye_s_m(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end

for j = 2:length(xe2_s)-2
    for k = 1:length(ye_s_r)
        rectangle('Position',[xe2_s(j),ye_s_r(k),cd,cwe]);
        rectangle('Position',[xe2_s(j)-sd,ye_s_r(k),sd,cwe],'LineWidth',2,'LineStyle','--');
    end
end
for k = 1:length(ye_s_r2)
    rectangle('Position',[xe2_s(length(xe2_s)-1),ye_s_r2(k),cd,cwe]);
    rectangle('Position',[xe2_s(length(xe2_s)-1)-sd,ye_s_r2(k),sd,cwe],'LineWidth',2,'LineStyle','--');
end

%plotting exits
%rectangle('Position',[0,0,wall,ydim],'FaceColor','r')
rectangle('Position',[0,0,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[0,ydim-y_exit,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim-x_exit,0,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[xdim-x_exit,ydim-y_exit,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[x_exit+space_front+2*rwb*rb+space_back-x_exit/2,0,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[x_exit+space_front+2*rwb*rb+space_back-x_exit/2,ydim-y_exit,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[2*x_exit+2*space_front+2*rwb*rb+2*rwe*re1+2*space_back-x_exit/2,0,x_exit,y_exit],'FaceColor','g')
rectangle('Position',[2*x_exit+2*space_front+2*rwb*rb+2*rwe*re1+2*space_back-x_exit/2,ydim-y_exit,x_exit,y_exit],'FaceColor','g')
axis([0 xdim 0 ydim]);
%grid off
set(gcf, 'Position', [100 100 2.5*(xdim) 5*(ydim)]);
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

        
        swarm(i, 1, 1) = swarm(i, 1, 1) + 6.5*swarm(i, 2, 1)*dt;     %update x position
        swarm(i, 1, 2) = swarm(i, 1, 2) + 6.5*swarm(i, 2, 2)*dt;     %update y position
            
        x = swarm(i, 1, 1);
        y = swarm(i, 1, 2);
        
        for j = (i+1) : swarm_size
                dx = swarm(j,1,1)-swarm(i,1,1);
                dy = swarm(j,1,2)-swarm(i,1,2);
                
                distance = sqrt(dx*dx + dy*dy);
                
                if( distance < tol)
                
                    x2 = round((x+dx/2)*numx);
                    y2 = round((y+dy/2)*numy);
                    %chair_repel = pc(x2,y2);
                
                    u = c3*(Ca/La * exp(-distance / La) - Cr/Lr * exp(-distance / Lr));
                    %u = -c3/r^2/chair_repel;
                    fx(i) = fx(i) + u*dx/distance;
                    fy(i) = fy(i) + u*dy/distance;
                    fx(j) = fx(j) - u*dx/distance;
                    fy(j) = fy(j) - u*dy/distance; 
                end
         end
 
%         f(i) = sqrt(fx(i)*fx(i)+fy(i)*fy(i));
%         f_count = f_count+1;
%         [ f_count,f_max,f_ave,bin_f ] = f_data( f(i), f_count, f_sum,f_max,partition_f, bin_f );
        
        
        ix=round(x*numx);
        iy=round(y*numy);
        %swarm(i,4,1)=p(ix,iy); % fitness evaluation

        
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
          % the order of exits get randomized every iter, so that agents
          % would choose a different exit when min(dis) to differnt exit is
          % equal
          exits = exits(randperm(8),:);
          for j=1:8  
            xdis(j) = exits(j,1)-swarm(i,1,1);
            ydis(j) = exits(j,2)-swarm(i,1,2);
            dis(j) = xdis(j)*xdis(j) + ydis(j)*ydis(j);
            
          end
          [r2,pos]=min(dis);
          r = sqrt(r2);
         
          
          % if at the exit, add index to a vector called deletions
         
          if r<=7.5*rad
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
%         if swarm_size <= 0.05*291
%             total_iter = iter;
%         end
        
    
    
    % Plotting the swarm
    
    if mod(iter,20)==0
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
         set(gcf, 'Position', [100 100 2.5*(xdim) 5*(ydim)]);
         frame = getframe;
         %writeVideo(writerObj,frame);
        
        %pause(4)
    end
    
    if(swarm_size == 0) 
                break; end
end

hold off




%close(writerObj);
%% taken out the video output part
 