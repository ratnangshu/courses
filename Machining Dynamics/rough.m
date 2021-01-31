clc;
clear;

%% Given Data
Y = 750; % Yield Strength
alpha_n = 20*(pi/180); % normal rake angle
%alpha_n = 20*(pi/180); % normal rake angle
i = 40*(pi/180); % angle of obliquity
mu = 0.5; % coefficient of friction
h = 1; % uncut chip thickness
hc = 2.2; % cut chip thickness
b = 10; % width of cut

%% 
%beta_a = atan(mu); % friction angle
beta_a = 34.6*(pi/180);

% syms phi_n phi_i eta theta_i theta_n 
% 
% eqn1 = sin(theta_i) == sin(beta_a)*sin(eta);
% eqn2 = tan(theta_n + alpha_n) == tan(beta_a)*cos(eta);
% eqn3 = tan(eta) == (tan(i)*cos(phi_n - alpha_n) - cos(alpha_n)*tan(phi_i))/(sin(phi_n));
% 
% %% Minimum Power
% Pt_p = (cos(theta_n) + tan(theta_i)*tan(i))/((cos(theta_n + phi_n)*cos(phi_i) + tan(theta_i)*sin(phi_i))*sin(phi_n));
% eqn4 = diff(Pt_p, phi_n) == 0;
% eqn5 = diff(Pt_p, phi_i) == 0;
v = 0.1;
for ang = 1:2:40
%for ang = 35
    i = ang*pi/180;
    eta_temp = i;
    for k = 0:10
        % Force Relation
        % SF = solve([eqn1, eqn2], [theta_i, theta_n]);
        eta_val = eta_temp;
        theta_i_val = asin(sin(beta_a)*sin(eta_val));
        theta_n_val_temp = atan(cos(eta_val)*tan(beta_a)) - alpha_n;
        theta_n_val = atan(sin(theta_n_val_temp)/cos(theta_n_val_temp));

        % Shear angle prediction
        syms phi_n phi_i
        Pt_p = (cos(theta_n_val) + tan(theta_i_val)*tan(i))/((cos(theta_n_val + phi_n)*cos(phi_i) + tan(theta_i_val)*sin(phi_i))*sin(phi_n));
        eqn4 = diff(Pt_p, phi_n) == 0;
        eqn5 = diff(Pt_p, phi_i) == 0;
        %SSA = vpasolve([eqn4, eqn5], [phi_i, phi_n], [-pi/2 pi/2; -pi/2 pi/2]);
        SSA = vpasolve([eqn4, eqn5], [phi_i, phi_n], [0; 40*pi/180]);
        phi_i_val = double(SSA.phi_i);
        phi_n_val = double(SSA.phi_n);

        % Velocity relation
        eta_new = atan(double(subs((tan(i)*cos(phi_n_val - alpha_n) - cos(alpha_n)*tan(phi_i_val))/(sin(phi_n_val)))));
        eta_temp = v*eta_temp + (1-v)*eta_new;
        
        % Power 
        Pt_p = (cos(theta_n_val) + tan(theta_i_val)*tan(i))/((cos(theta_n_val + phi_n_val)*cos(phi_i_val) + tan(theta_i_val)*sin(phi_i_val))*sin(phi_n_val));
        %plot(k,Pt_p,'bo'); hold on;
    end
    subplot(2,2,1)
    plot(ang,theta_n_val*(180/pi),'b.'); hold on;
    subplot(2,2,2)
    plot(ang,theta_i_val*(180/pi),'b.'); hold on;
    subplot(2,2,3)
    plot(ang,phi_n_val*(180/pi),'b.'); hold on;
    subplot(2,2,4)
    plot(ang,eta_temp*(180/pi),'b.'); hold on;
end

% %S = vpasolve([eqn1, eqn2, eqn3, eqn4, eqn5], [theta_i, theta_n, eta, phi_i, phi_n], [-pi pi; -pi pi; -pi pi; -pi pi; -pi pi]);
% S = solve([eqn1, eqn2, eqn3, eqn4, eqn5], [theta_i, theta_n, phi_i, phi_n]);
% 
% %eqn6 = eta == i; % Stabler's Rule
% S = vpasolve(subs([eqn1, eqn2, eqn4, eqn5, eqn6], eta, i), [theta_i, theta_n, phi_i, phi_n], [-pi pi; -pi pi; -pi pi; -pi pi]);
