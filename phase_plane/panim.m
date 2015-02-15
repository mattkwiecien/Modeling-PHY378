function panim
% x and y coordiantes made from unit circle
ppi = 0:.05:2*pi;
a = 2.*(cos(ppi))+.5;
b = 2.*(sin(ppi))+.5;

nframe = length(ppi);
mov(1:nframe) = struct('cdata',[],'colormap',[]);

fig = figure;
set(fig, 'Position', [250, 250, 1500, 500])
hold on

a_parabola = linspace(-3,1,40);
b_parabola = (a_parabola-1).^2 ./4 + a_parabola;

for i = 1:length(ppi)
    
    %Uses pando's code to grab the phase plane as well as a single solution
    %to the system of odes
    [x,y,px,py,xsol,ysol] = phplane(a(i),b(i));
    
    %plots the unit circle - most of this is just plotting formatting
    subplot(1,2,2)
    hold on
    grid on
    %plots the current values
    plot(a(i),b(i),'ro','markersize',8)
    %plots circle
    plot(a,b,'-r')
    line([-10 10], [-10 10], 'color', 'b', 'linewidth', 2)
    line([1 1], [-10 10], 'color', 'b', 'linewidth', 2)
    plot(a_parabola,b_parabola,'b', 'linewidth', 2)
    %creates a line to current values in green
    line([0 a(i)], [b(i) b(i)], 'color', 'g', 'linewidth', 2)
    line([a(i) a(i)], [0 b(i)], 'color', 'g', 'linewidth', 2)
    
    ylim([-2 3])
    xlim([-2 3])
    xL = xlim;
    yL = ylim;
    %x and y axis
    line(xL, [0 0], 'color', 'k', 'linewidth', 1)
    line([0 0],yL, 'color', 'k', 'linewidth', 1)
    
    %Creates a text box and puts it on the plot window based on what our
    %conditions are
    %
    % THIS IS WHERE WE HAVE TO FIGURE OUT THE CORRECT VALUES
    if b(i) > ((a(i)+4)^2 / 4) +a(i)
        legend('2 Complex Eigenvalues')
        if a(i)>1
            legend('2 Complex Eigenvalues | Spiral Source')
        elseif a(i)<1
            legend('2 Complex Eigenvalues | Spiral Sink')
        else
            legend('2 Complex Eigenvalues | Center')
        end    
    elseif b(i) < ((a(i)+4)^2 / 4)+a(i)
        legend('2 Real Eigenvalues')
        if a(i)>1
            legend('2 Real Eigenvalues | Source')
            if b(i)<1
                legend('2 Real Eigenvalues | Saddle')
            end
        elseif a(i)<1
            legend('2 Real Eigenvalues | Sink')
            if b(i)<1
                legend('2 Real Eigenvalues | Saddle')
            end
        end
    elseif b(i) == a(i)
        if a(i) > 1
            legend('2 Real Eigenvalues | Line of Source Points')
        elseif a(i) < 1
            legend('2 Real Eigenvalues | Line of Sink Points')
        end
    else 
        legend('Repeated Eigenvalues')
        if a(i)<1
            legend('Repeated Eigenvalues | Line Sink')
        elseif a(i)>1
            legend('Repeated Eigenvalues | Line Source')
        else
            legend('Repeated Eigenvalues | Center')
        end
    end
    
    title('The Unit Circle for values of a and b.')
    xlabel('a')
    ylabel('b')
    
    %plots the phase plane and the paarticular solution 
    subplot(1,2,1)
    hold on
    title('Phase Plane and Single Solution for system of ODEs at a,b.')
    plot(xsol,ysol,'m','linewidth',3)
    quiver(x,y,px,py)
    axis([-10 10 -10 10]);
    pause(.1)
    
    % ignore this - just for making Gifs
    % ---------------------------------------------------------------
    
    set(gcf,'color','w');
    mov(i) = getframe(gcf);
%     drawnow;
%     frame = getframe(1);
%     im=frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     outfile = 'test.gif';
%     
%     if i==1
%         imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
%     else
%         imwrite(imind,cm,outfile,'gif','WriteMode','append');
%     end
     % ---------------------------------------------------------------

    clf
    
end

movie2avi(mov,'plane.avi', 'compression', 'None');


end