function nonlinear(solplot)
%Plots phase plane, nullclines and some solutions
%   Uses lode.m

xnot=2;
ynot=2;
tmin= 0;
tmax = 100;
figure

[x,y] = meshgrid(-5:.35:5);
dx=  x + y - x.^3;
dy = -.5.*x;
r = (dx.^2 + dy.^2).^(0.5);
px = dx./r;
py = dy./r;

 quiver(x,y,px,py)

axis([0 3 0 3]);
 hold on
if solplot >0
[trk4,yrk4]=ode45(@lode,[tmin tmax],[xnot ynot]);
plot(yrk4(:,1),yrk4(:,2),'-k','LineWidth',2)
legend('Phase Plane','Solution for x_o = 2, y_o=2')
end
h=ezplot('x^3 - x',[0,3,0,3]);%,'-r','LineWidth',2);
set(h,'color','red','linewidth',2);
h=ezplot('0',[0,3,0,3]);%,'-r','LineWidth',2);
set(h,'color','magenta','linewidth',2)
end

