% data for business cabin
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

%get the plane potential function

[ p,numx,numy,xdim,ydim,x,y ]=Boeing737_plane(r,c,a,cw,rw,sd,cd,lr,wall,nc,as,rs,cs,ws,space_front,space_back,space_above,space_below,x_exit,y_exit);
[p_x,p_y,sz] = plane_grad( p,rad,numx,numy );


%plotting the plane function and grad
  figure;
  set(gcf, 'Position', [20 20 3000 900])
  mesh(y,x,p);
  axis vis3d
  pause(2)
  
  
  figure;
    spacing_x = 8;
    spacing_y = 6;
    starty = 5;
    quiver(y(starty:spacing_y:sz(2)),x(1:spacing_x:sz(1)),-p_y(1:spacing_x:sz(1),starty:spacing_y:sz(2)),-p_x(1:spacing_x:sz(1),starty:spacing_y:sz(2)));
    axis([0 ydim 0 xdim]);
    pause(2)