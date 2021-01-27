%% PHASE PLOT 
clc;
clear;

f = @(t,x) [-sin(x(1))*(1+x(2)); x(1) + x(2)];
%f = @(t,x) [x(2); x(1) - x(1)^2];

spacing = 0.4;
[X,Y] = meshgrid(-4:spacing:4);

U = zeros(size(X));
V = zeros(size(X));

t=0; 
for i = 1:numel(X)
    xdot = f(t,[X(i); Y(i)]);
    U(i) = xdot(1);
    V(i) = xdot(2);
end

quiver(X,Y,U,V,1); 
xlabel('y_1')
ylabel('y_2')
axis square;

hold on
[x10, x20] = ginput(1);
tim_end = 4;
[ts,xs] = ode45(f,[0,tim_end],[x10;x20]);
plot(xs(:,1),xs(:,2),'r')
plot(xs(1,1),xs(1,2),'bo') % starting point
plot(xs(end,1),xs(end,2),'ks') % ending point
xlim([-5 5])
ylim([-5 5])
hold off