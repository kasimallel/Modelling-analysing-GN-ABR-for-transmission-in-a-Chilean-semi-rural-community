%Graphs
%http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
%t = 1:1:150;
%t= t(:);
%Population size over time
figure(2)
plot(t,[(y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5)+y(:,6))], '-','LineWidth', 2.5, 'color', 'r')
hold on
plot(t,y(:,8)+y(:,7), '--', 'LineWidth', 2.5, 'color', 'r')
set(gca, 'FontSize', 16);
xlabel('Time')
ylabel('Population')
title(['Population size over time'] )
legend('Total population', 'Number of deaths attributed to BS')

%Infections by type of Bloodstream infection
figure(3)
plot(t,y(:,5), '-','LineWidth', 2.5, 'color', 'r')
hold on
plot(t,y(:,6), '--', 'LineWidth', 2.5, 'color', 'r')
set(gca, 'FontSize', 16);
xlabel('Time')
ylabel('Number of blooodstream infections')
title(['Bloodstream infections over time'] )
legend('Attributed to susceptible bacteria','Attributed to resistant bacteria')


figure(4)
plot(t,[(y(:,4)], '-','LineWidth', 2.5, 'color', 'r')
set(gca, 'FontSize', 16);
xlabel('Time')
ylabel('Population')
title(['Population size over time'] )
legend('Total population')



figure(4)
plot(t,R_n2, '-','LineWidth', 3.5, 'color', 'g')
set(gca, 'FontSize', 20);
xlabel('Time in days')
ylabel('R_0 values')
title([''] )
legend('R_0')




x=y(:,4);
y=y(:,6);
%t;
% interpolation on regular grid
xlin=linspace(min(x),max(x),50);
ylin=linspace(min(y),max(y),50);
[X,Y]=meshgrid(xlin,ylin);
Z=griddata(y,x,t,Y,X); 

% visualization
mesh(Z,Y,X);
axis tight; hold on
plot3(t,y,x,'.', 'MarkerSize',15)
xlabel('x')
ylabel('y')
title('Resistance and BS infections over time')
%surf(X,Y,Z)



figure(1)
%plot(y(:,6),y(:,4), 'LineWidth', 2.5, 'color', 'k')
p = polyfit(y(:,6),y(:,4),2);
plot(y(:,4),p, 'LineWidth', 2.5, 'color', 'k')
set(gca, 'FontSize', 11);
xlabel('Number of cases with BS infections caused by resistant bacteria')
ylabel('Number of individuals colonized by resistant bacteria')
title([''] )
legend('')

[population2,gof] = fit(y(:,6),y(:,4),'poly3');
plot(population2,y(:,6),y(:,4));
legend('Location','NorthWest');



