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

% repulsion strength
c3 = 1;%80.0;%40;35;%65


%get the plane potential function

[ p,numx,numy,xdim,ydim,x,y ]=Boeing777_plane_block2(rb,re1,re2,cb,ce,ab,ae,cwb,cwe,rwb,rwe,sd,cd,lrb,lre,ncb,nce,as,rs,cs,space_front,space_back,x_exit);
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