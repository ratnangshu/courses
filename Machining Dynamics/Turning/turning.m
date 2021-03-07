clc;
clear;

D = 57; % steel shaft diameter
alpha_f = 5*(pi/180); % side rake angle
alpha_p = -5*(pi/180); % back rake angle
approach_angle = 0; 
r = 0.8; % nose radius
a = 1; % radial depth of cut
c = 0.06; % feed rate
V = 240; % cutting speed

%% Region 1
h = c; % chip thickness is constant and equal to feed rate
psi_r = 0; % Side cutting edge angle,
% Angles of interest
i = atan(tan(alpha_p)*cos(psi_r) - tan(alpha_f)*sin(psi_r)); % Equivalent oblique angle,
alpha_o = atan(tan(alpha_f)*cos(psi_r) + tan(alpha_p)*sin(psi_r)); % Orthogonal rake angle,
alpha_n = atan(tan(alpha_o)*cos(i)); % Normal rake angle,


fprintf('Equivalent oblique angle in Region I: %d\n',i);
fprintf('Orthogonal rake angle in Region I: %d\n',alpha_o);
fprintf('Normal rake angle in Region I: %d\n\n',alpha_n);


% orthogonal parameters for the same steel and carbide tool with a rake angle alpha_r=−5°
tau_s = 1400*h + 0.327*V + 507;
beta_a = (33.69 - 12.16*h - 0.0022*V)*(pi/180);
rc = 2.71*h + 0.00045*V + 0.227;

% Cutting force coefficients using simplified armarego's model
eta = i; % simplified armarego model for analytical solution
beta_n = atan(tan(beta_a)*cos(eta));
phi_n = atan((rc*(cos(eta)/cos(i))*cos(alpha_n))/(1-rc*(cos(eta)/cos(i))*sin(alpha_n)));
eta = atan((cos(alpha_n)*tan(i))/(tan(phi_n + beta_n)) + (sin(alpha_n)*tan(i)));

% Force coefficients
Ktc = (tau_s/sin(phi_n))*((cos(beta_n-alpha_n)+tan(i)*tan(eta)*sin(beta_n))/(sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)^2)));
Kfc = (tau_s/(sin(phi_n)*cos(i)))*((sin(beta_n-alpha_n))/(sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)^2)));
Krc = (tau_s/sin(phi_n))*((cos(beta_n-alpha_n)*cos(beta_n-alpha_n)*tan(i)-tan(eta)*sin(beta_n))/(sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)^2)));

% Forces
Ft1 = Ktc*c*(a-r); Fx1 = Ft1;
Fr1 = Krc*c*(a-r); Fy1 = Fr1;
Ff1 = Kfc*c*(a-r); Fz1 = Ff1;

fprintf('F_x in Region I: %d N \n',Fx1);
fprintf('F_y in Region I: %d N \n',Fy1);
fprintf('F_z in Region I: %d N \n\n',Fz1);


%% Region 2
% distribution of chip thickness
theta0 = pi - acos(c/(2*r)); % total angular contact
theta = 0:0.01:theta0;
gamma = theta - asin((c/r)*sin(pi-theta));
h = r - sqrt(c^2 + r^2 - 2*c*r*cos(gamma)); % Instantaneous chip thickness at theta
figure 
plot(theta,h);
grid on;
xlabel('theta [degrees]');
ylabel('h [mm]');

psi_r = theta; % Side cutting edge angle,
% Angles of interest
i = atan(tan(alpha_p)*cos(psi_r) - tan(alpha_f)*sin(psi_r)); % Equivalent oblique angle,
alpha_o = atan(tan(alpha_f)*cos(psi_r) + tan(alpha_p)*sin(psi_r)); % Orthogonal rake angle,
alpha_n = atan(tan(alpha_o).*cos(i)); % Normal rake angle,

% orthogonal parameters for the same steel and carbide tool with a rake angle alpha_r=−5°
tau_s = 1400*h + 0.327*V + 507;
beta_a = (33.69 - 12.16*h - 0.0022*V)*(pi/180);
rc = 2.71*h + 0.00045*V + 0.227;

% Cutting force coefficients using simplified armarego's model
eta = i; % simplified armarego model for analytical solution
beta_n = atan(tan(beta_a).*cos(eta));
phi_n = atan((rc.*(cos(eta)./cos(i)).*cos(alpha_n))./(1-rc.*(cos(eta)./cos(i)).*sin(alpha_n)));
eta = atan((cos(alpha_n).*tan(i))./(tan(phi_n + beta_n)) + (sin(alpha_n).*tan(i)));

% Force coefficients
Ktc = (tau_s./sin(phi_n)).*((cos(beta_n-alpha_n)+tan(i).*tan(eta).*sin(beta_n))./(sqrt(cos(phi_n+beta_n-alpha_n).^2 + tan(eta).^2.*sin(beta_n).^2)));
Kfc = (tau_s./(sin(phi_n).*cos(i))).*((sin(beta_n-alpha_n))./(sqrt(cos(phi_n+beta_n-alpha_n).^2 + tan(eta).^2.*sin(beta_n).^2)));
Krc = (tau_s./sin(phi_n)).*((cos(beta_n-alpha_n).*cos(beta_n-alpha_n).*tan(i)-tan(eta).*sin(beta_n))./(sqrt(cos(phi_n+beta_n-alpha_n).^2 + tan(eta).^2.*sin(beta_n).^2)));

figure 
subplot(1,3,1)
plot(theta,Ktc)
xlim([0, theta0]);
grid on;
axis square
xlabel('theta [degrees]');
ylabel('K_{tc} [N/mm^2]');

subplot(1,3,2)
plot(theta,Kfc)
xlim([0, theta0])
grid on;
axis square
xlabel('theta [degrees]');
ylabel('K_{fc} [N/mm^2]');

subplot(1,3,3)
plot(theta,Krc)
xlim([0, theta0])
grid on;
axis square
xlabel('theta [degrees]');
ylabel('K_{rc} [N/mm^2]');

% Forces
Ft2 = 0; Fr2 = 0; Ff2 = 0;
Fx2 = 0; Fy2 = 0; Fz2 = 0;
dtheta = 0.01;
for th = 1:length(theta)
    dFt = Ktc(th)*h(th)*r*dtheta;
    dFr = Krc(th)*h(th)*r*dtheta;
    dFf = Kfc(th)*h(th)*r*dtheta;
    
    dFx = dFt;
    dFy = -dFf*sin(theta(th)) + dFr*cos(theta(th));
    dFz = dFf*cos(theta(th)) + dFr*sin(theta(th));
    
    % tangential, radial, and feed cutting forces
    Ft2 = Ft2 + dFt;
    Fr2 = Fr2 + dFr;
    Ff2 = Ff2 + dFf;
    
    % total cutting forces
    Fx2 = Fx2 + dFx;
    Fy2 = Fy2 + dFy;
    Fz2 = Fz2 + dFz;
end

fprintf('F_t in Region II: %d N \n',Ft2);
fprintf('F_r in Region II: %d N \n',Fr2);
fprintf('F_f in Region II: %d N \n\n',Ff2);

fprintf('F_x in Region II: %d N \n',Fx2);
fprintf('F_y in Region II: %d N \n',Fy2);
fprintf('F_z in Region II: %d N \n\n',Fz2);

Fx = Fx1 + Fx2;
Fy = Fy1 + Fy2;
Fz = Fz1 + Fz2;

fprintf('Total F_x: %d N \n',Fx);
fprintf('Total F_y: %d N \n',Fy);
fprintf('Total F_z: %d N \n\n',Fz);


% Torque
T = Fx*(D-a)/2;
fprintf('Torque T: %d Nmm \n\n',T);

% Power
P = Fx*V/60;
fprintf('Power P: %d W \n\n',P);


