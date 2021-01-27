clc;
clear;

syms x1 x2

f1 = -sin(x1)*(1+x2); % x_1_dot = f1
f2 = x1 + x2; % x_2_dot = f2

S = solve([f1==0, f2==0], [x1 x2]);
x_1 = S.x1;
x_2 = S.x2;

q_star = [x_1,x_2; pi,-pi]; % some equilibrium points

D = [diff(f1,x1), diff(f1,x2); diff(f2,x1), diff(f2,x2)]; % Jacobian matrix
for i=1:3
    [x1, x2] = deal(q_star(i,1), q_star(i,2));
    D_i(:,:,i) = double(subs(D)); % Calculating jacobian Matrix at each of the equilibrium points
end

i=3;
[v,d] = eig(D_i(:,:,i)); % Eigenvalues and Eigenvectors at those points

%% PHASE PLOT 

f = @(t,x) [-sin(x(1))*(1+x(2)); x(1) + x(2)];

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

quiver(X,Y,U,V,1); % quiver directions showing the direction of phase plot
xlabel('x_1')
ylabel('x_2')
axis square;
xlim([-4 4]);
ylim([-4 4]);

hold on

% some standard lines showing the phase plot
for x10 = -4:0.8:4
    for x20 = -4:0.8:4
        tim_end = 4;
        [ts,xs] = ode45(f,[0,tim_end],[x10;x20]);
        plot(xs(:,1),xs(:,2),'r');
    end
end

% some lines showing the phase plot in the region near the fixed points
e = 0.1;
for xq = [-pi, 0, 1, pi]
    for x10 = -e:e/2:e
        for x20 = -e:e/2:e
            tim_end = 10;
            [ts,xs] = ode45(f,[0,tim_end],[x10+xq;x20-xq]);
            plot(xs(:,1),xs(:,2),'r');
        end
    end
end



