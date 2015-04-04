% rerun 777 to get data for Boeing 777


    % Initialize random number generator
random_number_generator = rng('shuffle');
rand(1,100);
v_max=0;

% Parameters
    iterations = 4000;
    w = 0.85; %inertia factor
    c1 = 1.5;%2.5; %self-confidence factor
    c2 = 0.2; %swarm confidence

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
    [ p,numx,numy,xdim,ydim ]=Boeing777_plane_block1(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit);
    %adding chair
    [ pc]=chair_Boeing777(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit);
    %gradiant for particle local search
    [p_x,p_y] = plane_grad( p,rad,numx,numy );
    %[pu_x,pu_y] = plane_grad_unscaled( p,rad,numx,numy );

    
    
    % Atttraction length
    La = 48;%36
    % Repultion length
    Lr = 24;%18
    % Attraction strength
    Ca = 5; %30
    % Repultion strength
    Cr = 15;%50

    
    % Here make a matrix of all global best positions (8 in this case)
    % call that matrix exits and then have each row be an x,y position of
    % an exit
    exits = zeros(8,2);
%     exits(1,1) = x_exit;
%     exits(2,1) = x_exit;
    exits(3,1) = xdim-x_exit;
    exits(4,1) = xdim-x_exit;
    exits(5,1) = x_exit+space_front+2*rb*rwb+space_back+0.5*x_exit;
    exits(6,1) = x_exit+space_front+2*rb*rwb+space_back+0.5*x_exit;
    exits(7,1) = 2.5*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back;
    exits(8,1) = 2.5*x_exit+2*space_front+2*rb*rwb+2*re1*rwe+2*space_back;
%     exits(1,2) = y_exit/2;
%     exits(2,2) = ydim-y_exit/2;
    exits(3,2) = y_exit/2;
    exits(4,2) = ydim-y_exit/2;
    exits(5,2) = y_exit/2;
    exits(6,2) = ydim-y_exit/2;
    exits(7,2) = y_exit/2;
    exits(8,2) = ydim-y_exit/2;


    total_iter = []; 
for rerun = 1:5
    % Initialization
    
    swarm_size = 19+158+114;

    %preset the matrix
    swarm = zeros(swarm_size, 4, 2);
    s = zeros(swarm_size, 4, 2);


    
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

                    r = sqrt(dx*dx + dy*dy);

                    if( r < tol)

                        x2 = round((x+dx/2)*numx);
                        y2 = round((y+dy/2)*numy);
                        %chair_repel = pc(x2,y2);

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
              
              if v(i)>=v_max
                v_max = v(i);
              end

        end
            %deleting the agent who exits the plane
            swarm(deletions,:,:)=[];
            swarm_size=length(swarm(:,1,1));
            if swarm_size <= 0.0*291
                total_iter(end+1) = iter;
            end

            if(swarm_size <= 0.0*291) 
                break; end

    end

end

ave_iter = sum(total_iter)/length(total_iter)
v_max
% close(writerObj);
%% taken out the video output part
 