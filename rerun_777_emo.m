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
total_iter = [];
for rerun = 1:3
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
    %width of the exits
    x_exit = 80;%20;
    %length of the exits
    y_exit = 10;
    %for deleting agents
    deletions = [];
    %thickness of the wall
    wall=0;
    wall_moving=zeros(1,iterations);
    wall_moving(1) = y_exit;
    fire_rate = 0.0;
    fire_accel = 0.0000005;
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
    [ p,numx,numy,xdim,ydim ]=Boeing777_plane_block4(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit);
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

    %%========================March 17th edits===============================%%
    %for random selecting a fraction of agents to be panic and freeze
    fraction = 0.1;%0.1;%0.2;
    rand_indices = randperm(swarm_size);
    swarm(rand_indices(1:floor(fraction*swarm_size)),3,2)=ones(fraction*swarm_size,1);  
    % this can be random later like 0.5 + rand(1:floor(fraction*swarm_size))*0.5;

    % %hardcode for visual effect
    % swarm(85,3,2) = rand;

    %setting parameter for panic susceptibility
    ps = 0.5;
    emo_thre = 0.7; %emotion threshold (when panic level>emo_thre*swarm(i,3,1)->freeze)

    
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



    %%outputting video
    %  writerObj = VideoWriter('plane.avi');
    %  writerObj = VideoWriter('C:\Users\Junyuan Lin\Dropbox\Pepp\plane project\plane_side', 'MPEG-4');
    % % %witerObj = VideoWriter('plane.avi','Uncompressed AVI');
    %  writerObj.FrameRate = 4;
    %  writerObj.Quality= 100;
    %  open(writerObj)


    % Here make a matrix of all global best positions (8 in this case)
    % call that matrix exits and then have each row be an x,y position of
    % an exit
    exits = zeros(8,2);
    exits(1,1) = x_exit/2.0;
    exits(2,1) = x_exit/2.0;
%     exits(3,1) = xdim-x_exit/2.0;
%     exits(4,1) = xdim-x_exit/2.0;
    exits(5,1) = x_exit+space_front+2*rb*rwb+space_back+0.5*x_exit;
    exits(6,1) = x_exit+space_front+2*rb*rwb+space_back+0.5*x_exit;
    exits(7,1) = 2.5*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back;
    exits(8,1) = 2.5*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back;
    exits(1,2) = y_exit/2;
    exits(2,2) = ydim-y_exit/2;
%     exits(3,2) = y_exit/2;
%     exits(4,2) = ydim-y_exit/2;
    exits(5,2) = y_exit/2;
    exits(6,2) = ydim-y_exit/2;
    exits(7,2) = y_exit/2;
    exits(8,2) = ydim-y_exit/2;






    % Iterations
    iter = 1;
    for iter = 1: iterations
        fx = zeros(swarm_size,1);
        fy = zeros(swarm_size,1); 
        fq = zeros(swarm_size,1);

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

            swarm(i, 4, 2) = swarm(i, 3, 2)*swarm(i, 3, 1) > emo_thre;
                % see whether the panic level is above the threshold
                % swarm(i,4,2) is a boolean variable
            
            for j = (i+1) : swarm_size
                    dx = swarm(j,1,1)-swarm(i,1,1);
                    dy = swarm(j,1,2)-swarm(i,1,2);

                    distance = sqrt(dx*dx + dy*dy);

                    dq = 0;
                    dq_f = 0; %emotion difference for updating force between agents

%=============================revised dq update==========================%
                if distance < 1.5*tol
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
                

                        x2 = round((x+dx/2)*numx);
                        y2 = round((y+dy/2)*numy);
                        %chair_repel = pc(x2,y2);

                        u = c3*(Ca/La * exp(-distance / La) - Cr/Lr * exp(-distance / Lr));
                        %u = -c3/distance^2/chair_repel;
                        fx(i) = fx(i) + u*dx/distance/(1-dq_f);
                        fy(i) = fy(i) + u*dy/distance/(1-dq_f);
                        fx(j) = fx(j) - u*dx/distance/(1+dq_f);
                        fy(j) = fy(j) - u*dy/distance/(1+dq_f); 
                end
%========================================================================%               
              
                    
                    if abs(sqrt((swarm(i,1,1)-(xdim-x_exit/2.0))*(swarm(i,1,1)-(xdim-x_exit/2.0))+(swarm(i,1,2)-(ydim-wall_moving(iter)))*(swarm(i,1,2)-(ydim-wall_moving(iter)))))<15*rad
                        dq=0.95-swarm(i,3,1);
                        fq(i) = fq(i)+0.5*dq;
                    end
                    if abs(sqrt((swarm(i,1,1)-(xdim-x_exit/2.0))*(swarm(i,1,1)-(xdim-x_exit/2.0))+(swarm(i,1,2)-(wall_moving(iter)))*(swarm(i,1,2)-(wall_moving(iter)))))<15*rad
                        dq=0.95-swarm(i,3,1);
                        fq(i) = fq(i)+0.5*dq;
                    end
                    %2.5*x_exit+2*rb*rwb+2*space_front+2*space_back+2*re1*rwe
                    %1.5*x_exit+2*rb*rwb+space_front+space_back
            end
            wall_moving(iter+1)=wall_moving(iter)+fire_rate;
            fire_rate= fire_rate+fire_accel;
    %         f(i) = sqrt(fx(i)*fx(i)+fy(i)*fy(i));
    %         f_count = f_count+1;
    %         [ f_count,f_max,f_ave,bin_f ] = f_data( f(i), f_count, f_sum,f_max,partition_f, bin_f );


            ix=round(x*numx);
            iy=round(y*numy);
            %swarm(i,4,1)=p(ix,iy); % fitness evaluation


        end

        %deleting the agent who exits the plane
            swarm(deletions,:,:)=[];
            swarm_size=length(swarm(:,1,1));
            if swarm_size <= 0.0*291
                total_iter(end+1) = iter;
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
              r_agent = sqrt(r2);


              % if at the exit, add index to a vector called deletions

              if r_agent<=7.5*rad
                  deletions = i;
              end

    %           r_count = r_count+1;
    %           r_sum = r_sum+r_agent;
              %[ r_max,r_ave ] = distance_data( r_agent,r_max,r_count,r_sum );



              % This is more efficient than gbest function.  
              gbest_x = xdis(pos)/(100.0+r_agent);
              gbest_y = ydis(pos)/(100.0+r_agent);

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
               v(i) = sqrt(swarm(i,2,1)*swarm(i,2,1)+swarm(i,2,2)*swarm(i,2,2));
                if v(i)>=v_max
                    v_max = v(i);
                end
               if fq(i) ~= 0
                     swarm(i, 3, 1)= swarm(i, 3, 1)+fq(i)/8*dt;  %8 is the nearest number
               else
                     swarm(i,3,1) = 0.999*swarm(i, 3, 1);
               end
             %for setting data for scatter plotting
    %             scatter_plot = zeros (swarm_size,3);
    %             scatter_plot(:,1) = swarm(:,3,1);         
    %             scatter_plot(:,3) = ones(swarm_size,1); 
    %             
        end


            
        if(swarm_size <= 0.0*291) 
                    break; end
    end
end

ave_iter = sum(total_iter)/length(total_iter)
v_max



% close(writerObj);
%% taken out the video output part
 