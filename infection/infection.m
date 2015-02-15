function infection(days)
%Parameters of Population

a = 3e-5;       % chance of infection when H and V interact
C = 4.6e-2;     % death rate of virus cells
bh = 0.515;     %Birth Rate of healthy cells 
bv = 5;         % birth rate of virus cells
tau = 4;        %latency period with virus but not infected

di = 4;         %death rate of infected cells
dn = 0.085;     % ratio between dn/di should be ~ .021

h0 = 1e6;       % healthy cells
i0 = 0;         % infected cells
v0 = 1;         % virus cells

% ========================================================================
% No Latency Period
% ========================================================================

tspan = linspace(0,days,1000);
xinit = [h0;i0;v0];

[T,X] = ode45(@dt,tspan,xinit);

H = X(:,1);
I = X(:,2);
V = X(:,3);

figure
hold on
set(gca,'fontsize',16)
plot(T,H,'g','linewidth', 3)
plot(T,I,'b','linewidth', 3)
plot(T,V,'r','linewidth', 3)
title('Model of Influenza A Spread for 50 Virus Cells in 200','fontsize',18)
xlabel('Days','fontsize',18)
ylabel('Number of Cells','fontsize',18)
legend('Healthy Cells', 'Infected Cells', 'Virus Cells')

% ========================================================================
% Latency Period
% ========================================================================
h1 = 1e6; % healthy cells
i1 = 0; % infected cells
v1 = 1; % virus cells
l1 = 0; % latent cells

xinit2 = [h1;i1;v1;l1];

[T2,X2] = ode45(@dt2,tspan,xinit2);

H2 = X2(:,1);
I2 = X2(:,2);
V2 = X2(:,3);
L2 = X2(:,4);

figure
hold on
set(gca,'fontsize',16)
plot(T,H2,'g','linewidth', 3)
plot(T,I2,'b','linewidth', 3)
plot(T,V2,'r','linewidth', 3)
plot(T,L2,'m','linewidth', 3)
title('Model of Influenza A Spread for 50 Virus Cells in 200 with Delay Period','fontsize',18)
xlabel('Days','fontsize',18)
ylabel('Number of Cells','fontsize',18)
legend('Healthy Cells', 'Infected Cells', 'Virus Cells', 'Virus Carrying but not Infected')

    function ddt = dt(t,x)
        hi = x(1);
        ii = x(2);
        vi = x(3);
        ddt(1) = -a*hi*vi + bh*(1 - (hi/h0))*hi;
        ddt(2) =  a*hi*vi - di*ii;
        ddt(3) =  bv*ii*(1-dn) - C*vi;
        if hi <=1
            ddt(1) = 0;
        end
        ddt = ddt';
    end
    function ddt2 = dt2(t2,x2)
        hi = x2(1);
        ii = x2(2);
        vi = x2(3);
        li = x2(4);
        ddt2(1) = -a*hi*vi + bh*(1 - (hi/h1))*hi;
        ddt2(2) =  tau*li - di*ii;
        ddt2(3) =  bv*ii*(1-dn) - C*vi;
        ddt2(4) =  a*hi*vi - tau*li;
        
        if hi <=1
            ddt2(1) = 0;
        end
        ddt2 = ddt2';
    end

end

