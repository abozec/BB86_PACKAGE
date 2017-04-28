    function draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
% ----------------------------------
% draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
%
% draws arrow with filled "head"
% Note: first X-coordinates x1, x2, then Y-coordinates, Y1,y2, ...
% Draws an arrow between pnt (x1,y1) and pnt (x2,y2)
% cf - scaling coefficient of the arrowhead (0 to 1)
% beta - angle between vector and arrow head beams (degrees)
% v_col - color ([R G B]) % default color is black
% lwd - line width, default = 1
% 
 if (isempty(v_col)); v_col=[0,0,0]; end; % default color is black
 if (isempty(lwd)); lwd=1.; end; % default line width

  hold on
  uu=x2-x1;
  vv=y2-y1;
  sp=sqrt(uu.*uu+vv.*vv);	  
  alfa=atan2(uu,vv);		      % vector angle from Y
  beta=beta*pi/180;
  var=cf*sp;				          % scaling of the arrow head
  dX2=var.*sin(alfa-beta);		  % arrow head coordinates
  dX3=var.*sin(alfa+beta);		  % 
  dY2=var.*cos(alfa-beta);
  dY3=var.*cos(alfa+beta);
  ax2=x2-dX2;
  ax3=x2-dX3;
  ay2=y2-dY2;
  ay3=y2-dY3;
%keyboard
  x2v=x1+(1-cf)*uu;
  y2v=y1+(1-cf)*vv;
  p1=plot([x1 x2v],[y1 y2v],'Color',v_col);     %vector
%  p2=plot([x2 ax2],[y2 ay2],'Color',v_col);		% arrow head
%  p3=plot([x2 ax3],[y2 ay3],'Color',v_col);
 
  set(p1,'linewidth',lwd);
%  set([p2,p3],'linewidth',0.5);

    X=[x2,ax2,ax3,x2];
    Y=[y2,ay2,ay3,y2];
    H=fill(X,Y,v_col);
    set(H,'edgecolor',v_col);


%  hold off
