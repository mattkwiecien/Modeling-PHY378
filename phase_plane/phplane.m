function [ x,y,px,py,xsol,ysol ] = phplane(a,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes her

%modified pando's code so it also produces a solution to our ode for a
%given initial coniditon
[x,y] = meshgrid(-10:1:10);

dx =  a.*x + b.*y;
dy = -1.*x - 1.*y;

r = (dx.^2 + dy.^2).^(0.5);

px = dx./r;
py = dy./r;

tspan = linspace(0,10,50);

[T,sol] = ode45(@ddt,tspan,[5; 5]);

xsol = sol(:,1);
ysol = sol(:,2);

    function dt = ddt(t,xx)
       
        dt(1) =  a*xx(1) + b*xx(2);
        dt(2) = -1*xx(1) - 1*xx(2);
        
        dt = dt';
        
    end

end

