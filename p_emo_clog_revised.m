% Particle Swarm Optimization Asiana Airplane Simulation 
% Simulates the movements of a swarm to minimize any functions with a 
% moving global minima
%  
% The swarm matrix is  
%
% swarm(index, [1-location, 2-velocity, 3-emotion, 4-best pos/panic level], 
%[1-x or the value component/agitated level(speeding up), 
%2-y component/panic susceptibility])

% Initialization
% Parameters
iterations = 2000;
w = 0.85; %inertia factor
c1 = 2.0; %self-confidence factor
c2 = 1; %swarm confidence
swarm_size = 60;


dt = 0.1; %time-step
% Number of rows
r = 10;
% Number of columns
c = 6; % 3 chairs, (1 for aisle), 3 chairs
% Aisle width, one chair width
a = 8;


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
wall=zeros(1,2000);
wall(1) = 20;
fire_rate = 0.0;
fire_accel = 0.0000005;
%space in the front of the plane
space_front = 21;
%space in the back of the plane
space_back = 250;
%person radius
rad = 9;

% aisle stength
as = 20;
% chair strength
cs = 20;
% row stength 
rs = 10;
%height of the wall
ws = 30;

% repulsion strength
repel = 15;%35;

% number of chairs per row
nc = 3;

[ p,numx,numy,xdim,ydim ]=plane(r,c,a,cw,rw,sd,cd,lr,wall(1),nc,as,rs,cs,ws,space_front,space_back);
[ pc]=chair_only(r,c,a,cw,rw,sd,cd,lr,wall(1),nc,as,rs,cs,ws,space_front,space_back);
[p_x,p_y] = plane_grad_unscaled( p,rad,numx,numy ); %here, since we've rescaled
%the gradient, so we are going to use plane_grad_unscaled funtion
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
swarm(:, 4, 2) = 0;             % initial panic level
swarm(:, 2, :) = 0;             % initial velocity
swarm(:, 3, :) = 0;             % initial emotion


%%========================March 6rd edits===============================%%
%for random selecting a fraction of agents to be panic and freeze
%fraction = 0.05;%0.1;%0.2;
% rand_indices = randperm(swarm_size);
%swarm(rand_indices(1:floor(fraction*swarm_size)),3,2)=ones(fraction*swarm_size,1);  % this can be random later like 0.5 + rand(1:floor(fraction*swarm_size))*0.5;


%hard code for visual effects
swarm(9, 3, 2) = 1;
swarm(21, 3, 2) = 1;
swarm(16, 3, 2) = 1;
%setting parameter for panic susceptibility
ps = 0.5;
emo_thre = 0.8; %emotion threshold (when panic level>emo_thre*swarm(i,3,1)->freeze)



% Atttraction length
La = 40;%36
% Repultion length
Lr = 20;%18
% Attraction strength
Ca = 20; %30
% Repultion strength
Cr = 50;%50

%% stuff for scaling force
% maxfx = zeros(1,iterations);
% maxfy = maxfx;
% minfx = maxfx;
% minfy = maxfx;
% Iterations

%% For outputting video
%  writerObj = VideoWriter('plane.avi');
%  writerObj = VideoWriter('C:\Users\Junyuan Lin\Desktop\Plane Project-clog\plane_emo_clog_3.25', 'MPEG-4');
% % %witerObj = VideoWriter('plane_clog_3.25.avi','Uncompressed AVI');
%  writerObj.FrameRate = 5;
%  writerObj.Quality= 100;
%  open(writerObj)


%% main loop
iter = 1;
while iter<iterations
    %min(swarm(:, 1, 1))<(xdim-190)
    fx = zeros(swarm_size,1);
    fy = zeros(swarm_size,1); 
    fq = zeros(swarm_size,1);
  
    for i = 1 : swarm_size
        
        s(i, 1, 1) = swarm(i, 1, 1) + swarm(i, 2, 1)*dt;     %update x position
        s(i, 1, 2) = swarm(i, 1, 2) + swarm(i, 2, 2)*dt;     %update y position
        
        if  s(i, 1, 1)>=xdim-rad      
            swarm(i, 2, 1) = 0;           
            swarm(i, 2, 2) = 0;
        end
        
        swarm(i, 1, 1) = swarm(i, 1, 1) + swarm(i, 2, 1)*dt;     %update x position
        swarm(i, 1, 2) = swarm(i, 1, 2) + swarm(i, 2, 2)*dt;     %update y position
              
        
        swarm(i, 4, 2) = swarm(i, 3, 2)*swarm(i, 3, 1) > emo_thre;
                % see whether the panic level is above the threshold
                % swarm(i,4,2) is a boolean variable
                
                
        x = swarm(i, 1, 1);
        y = swarm(i, 1, 2);
        
        for j = (i+1) : swarm_size
                
                dx = swarm(j,1,1)-swarm(i,1,1);
                dy = swarm(j,1,2)-swarm(i,1,2);
                
                r = sqrt(dx*dx + dy*dy);
                dq = 0;  %difference of emotion for updating swarm(i,3,1)
                dq_f = 0; %emotion difference for updating force between agents
                
                if r<1.5*rad
                    dq = swarm(j,3,1)-swarm(i,3,1);
                    dq_f = swarm(j,3,1)-swarm(i,3,1)*(1-swarm(i,4,2));
                    if  dq>0
                        fq(i) = fq(i)+6*dq;
                        fq(j) = fq(j)-3*dq;
                        dq_f = min(dq_f,0.5);
                    else
                        fq(i) = fq(i)+3*dq;
                        fq(j) = fq(j)-6*dq;
                        dq_f = max(dq_f,-0.5);
                    end     
                end
                
                x2 = round((x+dx/2)*numx);
                y2 = round((y+dy/2)*numy);
                chair_repel = pc(x2,y2);
                
                u = -repel/r^2/chair_repel;
                
                
                fx(i) = fx(i) + u*dx/r/(1-dq_f);
                fy(i) = fy(i) + u*dy/r/(1-dq_f);
                fx(j) = fx(j) - u*dx/r/(1+dq_f);
                fy(j) = fy(j) - u*dy/r/(1+dq_f);
                
               
                   
            if swarm(i,1,1)-wall(iter)<5*rad
                dq=1-swarm(i,3,1);
                fq(i) = fq(i)+0.5*dq;
            end
            
        end
        %fx(i)=fx(i)/swarm_size;
        %fy(i)=fy(i)/swarm_size;
        %fq(i)=fq(i)/swarm_size;
        wall(iter+1)=wall(iter)+fire_rate;
        fire_rate= fire_rate+fire_accel;
        
        ix=round(x*numx);
        iy=round(y*numy);
        swarm(i,4,1)=p(ix,iy); % fitness evaluation
        
%         if val <= swarm(i, 4, 1)                 % if new position is better
%             swarm(i, 3, 1) = swarm(i, 1, 1);    % update best x,
%             swarm(i, 3, 2) = swarm(i, 1, 2);    % best y postions
%             swarm(i, 4, 1) = val;               % and best value
%         end
    end

    [temp, gbestpos] = min(swarm(:, 4, 1));        % global best position
    
    %ux1 = rand();
    %ux2 = rand();
    %uy1 = rand();
    %uy2 = rand();
    % updating velocity vectors
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
          
          
          %%========the effects of panic factor===============%%
          if swarm(i, 4, 2)
              %reduce the overall speed of agents and the effects of the force exerted on them are more observable.
              swarm(i, 2, 1) = fx(i) + (w*swarm(i,2,1) - c1*ux1*p_x(ix,iy) + c2*ux2*gbest_x*(1+swarm(i,3,1)))*ps;
              swarm(i, 2, 2) = fy(i) + (w*swarm(i,2,2) - c1*uy1*p_y(ix,iy) + c2*uy2*gbest_y*(1+swarm(i,3,1)))*ps;
          else
              swarm(i, 2, 1) = fx(i) + w*swarm(i,2,1) - c1*ux1*p_x(ix,iy) + c2*ux2*gbest_x*(1+swarm(i,3,1));   
              swarm(i, 2, 2) = fy(i) + w*swarm(i,2,2) - c1*uy1*p_y(ix,iy) + c2*uy2*gbest_y*(1+swarm(i,3,1));
          end

          if fq(i) ~= 0
                swarm(i, 3, 1)= swarm(i, 3, 1)+fq(i)/8*dt;  %8 is the nearest number
          else
             swarm(i,3,1) = 0.999*swarm(i, 3, 1);
          end
          
          
          
    end
    
   % gbest(max(gx),max(gy))
    
    % Plotting the swarm
    if mod(iter,5)==1
        clf
        hold on
    %plot(swarm(:, 1, 1), swarm(:, 1, 2), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','MarkerSize',3*rad); % drawing swarm movements
%         for i = 1: swarm_size
%             scatter(swarm(i, 1, 1),swarm(i, 1, 2),500,[swarm(i, 3, 1) 0 1],'filled');
%         end
        
        %plotting the panic passengers -March 3rd (the panic ones will be cyan)
        for i = 1: swarm_size
            if ~swarm(i, 4, 2)
                scatter(swarm(i, 1, 1),swarm(i, 1, 2),500, [swarm(i, 3, 1) 0 1],'filled');
            else 
                scatter(swarm(i, 1, 1),swarm(i, 1, 2),500, [0 swarm(i, 3, 1) 1],'filled');
            end
        end
            
        quiver(swarm(:,1,1),swarm(:,1,2),swarm(:,2,1),swarm(:,2,2));
%     plot(10*pi,10*pi,'O','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     plot(30*pi,10*pi,'gO');
%     plot(22,32,'O','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     plot(108,58,'gO');
        for j = 1:length(x_s)
            for k = 1:length(y_s)
                rectangle('Position',[x_s(j),y_s(k),cd,nc*cw]);
            end
        end
        rectangle('Position',[0,0,wall(iter),ydim],'FaceColor','r')
        rectangle('Position',[xdim-210,nc*cw,20,a],'FaceColor','g')
        axis([0 xdim-190 0 ydim]);
        grid off
        set(gcf, 'Position', [100 100 2.5*(xdim-190) 5*ydim]);
        frame = getframe;
        %writeVideo(writerObj,frame);
         %pause(4)
    end
iter=iter+1;
end
% sprintf('Number of iterations taken for full evacuation is %d',iter);
%close(writerObj);
 