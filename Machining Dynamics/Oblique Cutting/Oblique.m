clc;
clear;

%% Given Data
Y = 750; % Yield Strength
alpha_n = 13*(pi/180); % normal rake angle
i = 8*(pi/180); % angle of obliquity
mu = 0.5; % coefficient of friction
h = 1; % uncut chip thickness
hc = 2.2; % cut chip thickness
b = 10; % width of cut

%% 
beta_a = atan(mu); % friction angle
eta_temp = i;

for k = 0:100
    syms phi_n phi_i eta theta_i theta_n 

    eqn1 = sin(theta_i) == sin(beta_a)*sin(eta);
    eqn2 = tan(theta_n + alpha_n) == tan(beta_a)*cos(eta);
    eqn3 = tan(eta) == (tan(i)*cos(phi_n - alpha_n) - cos(alpha_n)*tan(phi_i))/(sin(phi_n));

    %% Minimum Power
    Pt_p = (cos(theta_n) + tan(theta_i)*tan(i))/((cos(theta_n + phi_n)*cos(phi_i) + tan(theta_i)*sin(phi_i))*sin(phi_n));
    eqn4 = diff(Pt_p, phi_n) == 0;
    eqn5 = diff(Pt_p, phi_i) == 0;

    % Force Relation
    % SF = solve([eqn1, eqn2], [theta_i, theta_n]);
    syms eta
    theta_i = asin(sin(beta_a)*sin(eta));
    theta_n = atan(cos(eta)*tan(beta_a)) - alpha_n;

    eta = eta_temp;
    theta_i = double(subs(theta_i));
    theta_n = double(subs(theta_n));

    % Shear angle prediction
    SSA = vpasolve([subs(eqn4), subs(eqn5)], [phi_i, phi_n], [-pi pi; -pi pi]);
    phi_i = double(SSA.phi_i);
    phi_n = double(SSA.phi_n);

    % Velocity relation
    eta_temp = atan(double(subs((tan(i)*cos(phi_n - alpha_n) - cos(alpha_n)*tan(phi_i))/(sin(phi_n)))));

    % Power 
    Pt_p = (cos(theta_n) + tan(theta_i)*tan(i))/((cos(theta_n + phi_n)*cos(phi_i) + tan(theta_i)*sin(phi_i))*sin(phi_n));
    plot(k,Pt_p,'bo'); hold on;
end

% %S = vpasolve([eqn1, eqn2, eqn3, eqn4, eqn5], [theta_i, theta_n, eta, phi_i, phi_n], [-pi pi; -pi pi; -pi pi; -pi pi; -pi pi]);
% S = solve([eqn1, eqn2, eqn3, eqn4, eqn5], [theta_i, theta_n, phi_i, phi_n]);
% 
% %eqn6 = eta == i; % Stabler's Rule
% S = vpasolve(subs([eqn1, eqn2, eqn4, eqn5, eqn6], eta, i), [theta_i, theta_n, phi_i, phi_n], [-pi pi; -pi pi; -pi pi; -pi pi]);
