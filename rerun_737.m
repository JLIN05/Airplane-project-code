% run Boeing737_final multiple times to find average egress time

total_iter = [];

% Parameters
    iterations = 8000;
    w = 0.85; %inertia factor
    c1 = 1.5;%2.5; %self-confidence factor
    c2 = 0.5;%1.0; %swarm confidence

    % Atttraction length
    La = 48;%36
    % Repulsion length
    Lr = 24;%18
    % Attraction strength
    Ca = 5; %30
    % Repultion strength
    Cr = 15;%50





       

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

    

    dt = 0.05;%0.25; %time-step
    % Number of rows
    r = 7;
    % Number of columns
    c = 6; % 3 chairs, (1 for aisle), 3 chairs
    % Aisle width, one chair width
    a = 18;

    
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
    [ p,numx,numy,xdim,ydim ]=block_Boeing737_plane2(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back,space_above,space_below,x_exit,y_exit);
    %adding chair
    [ pc]=chair_Boeing737(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back,space_above,space_below,x_exit);
    %gradiant for particle local search
    [p_x,p_y] = plane_grad( p,rad,numx,numy );
    %[pu_x,pu_y] = plane_grad_unscaled( p,rad,numx,numy );
    
    
    
    % Here make a matrix of all global best positions (8 in this case)
    % call that matrix exits and then have each row be an x,y position of
    % an exit
    exits = zeros(8,2);
    exits(1,1) = space_out+wall+x_exit;
    exits(2,1) = space_out+wall+x_exit;
    exits(3,1) = xdim-space_out-wall-x_exit;
    exits(4,1) = xdim-space_out-wall-x_exit;
%     exits(5,1) = xdim/2.0-x_exit;
%     exits(6,1) = xdim/2.0-x_exit;
%     exits(7,1) = xdim/2.0+x_exit;
%     exits(8,1) = xdim/2.0+x_exit;
    exits(1,2) = space_below;%+y_exit/2;
    exits(2,2) = ydim-space_above;%-y_exit/2;
    exits(3,2) = space_below;%+y_exit/2;
    exits(4,2) = ydim-space_above;%-y_exit/2;
%     exits(5,2) = space_below;%+y_exit/2;
%     exits(6,2) = ydim-space_above;%-y_exit/2;
%     exits(7,2) = space_below;%+y_exit/2;
%     exits(8,2) = ydim-space_above;%-y_exit/2;
    
    
    
    %for figuring out the max velocity
    v_max = 0;
    
    % Iterations
for rerun = 1:5

       % Initialize random number generator
    random_number_generator = rng('shuffle');

    % Initialization
    
    swarm_size = 172;




    %preset the matrix
    swarm = zeros(swarm_size, 4, 2);
    s = zeros(swarm_size, 4, 2);


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



    dis = zeros(8,1);
    xdis = zeros(8,1);
    ydis = zeros(8,1);



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

                    distance = sqrt(dx*dx + dy*dy);

                    if( distance < tol)

                        x2 = round((x+dx/2)*numx);
                        y2 = round((y+dy/2)*numy);
                        chair_repel = pc(x2,y2);

                        u = c3*(Ca/La * exp(-distance / La) - Cr/Lr * exp(-distance / Lr));
                        %u = -c3/r^2/chair_repel;
                        fx(i) = fx(i) + u*dx/distance;
                        fy(i) = fy(i) + u*dy/distance;
                        fx(j) = fx(j) - u*dx/distance;
                        fy(j) = fy(j) - u*dy/distance; 
                    end
             end



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


              % here compute the distance to each global best.
              % then choose the gbest that is the shortest distance
              % then call that gbestpos(suggesting exit)
              for j=1:8  
                xdis(j) = exits(j,1)-swarm(i,1,1);
                ydis(j) = exits(j,2)-swarm(i,1,2);
                dis(j) = xdis(j)*xdis(j) + ydis(j)*ydis(j);
              end
              [r2,pos]=min(dis);
              r_agent = sqrt(r2);

              % if at the exit, add index to a vector called deletions

              if r_agent<=2.8*rad
                  deletions = i;
              end



              % This is more efficient than gbest function.  
              gbest_x = xdis(pos)/(100.0+r_agent);
              gbest_y = ydis(pos)/(100.0+r_agent);

              ux1 = rand();
              ux2 = rand();
              uy1 = rand();
              uy2 = rand();

              swarm(i, 2, 1) = fx(i) + w*swarm(i,2,1) - c1*ux1*p_x(ix,iy) + c2*ux2*gbest_x;   
              swarm(i, 2, 2) = fy(i) + w*swarm(i,2,2) - c1*uy1*p_y(ix,iy) + c2*uy2*gbest_y;
              v(i) = sqrt(swarm(i,2,1)*swarm(i,2,1)+swarm(i,2,2)*swarm(i,2,2));
              
              if v(i)>=v_max
                v_max = v(i);
              end
              
        end
            %deleting the agent who exits the plane
            swarm(deletions,:,:)=[];
            swarm_size=length(swarm(:,1,1));
            
            if swarm_size<= 0.0*172
                total_iter(end+1) = iter;
            end
            
            if(swarm_size <= 0.0*172) 
                break; end

       end
        
end



ave_iter = sum(total_iter)/length(total_iter)
v_max






% close(writerObj);
%% taken out the video output part
 