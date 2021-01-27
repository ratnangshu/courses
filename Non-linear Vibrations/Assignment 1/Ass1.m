clc;
clear;

diffeq1 = @(t,u)[u(2); -u(1)-u(1)^3];
diffeq2 = @(t,u)[u(2); -u(1)^3];

A = 0.001;
[t1,u1] = ode45(diffeq1,[0 10],[A; 0]);
[t2,u2] = ode45(diffeq2,[0 20000],[A; 0]);

%plot(t,u(:,1),'-o',t,u(:,2),'-o')
subplot(1,2,1)
plot(t1,u1(:,1),'-o')
xlabel('Time t');
ylabel('Solution 1');
axis square
grid on;

subplot(1,2,2)
plot(t2,u2(:,1),'-o')
xlabel('Time t');
ylabel('Solution 2');
axis square
grid on;



