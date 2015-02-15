function dy = lode(t,y)
%ODE system to be solved by odesystem
% y(1) = x; y(2) = y 

dy=zeros(2,1);
x = y(1);
yy = y(2);

dy(1) = x + yy - x^3;
dy(2) = -0.5*x;
% dy=exp(t)*sin(y);


end

