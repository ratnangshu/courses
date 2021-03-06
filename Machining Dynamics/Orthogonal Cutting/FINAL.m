clc;
clear;

%% Experimental Data
alpha_r = 0; % rake angle
b = 5; % width of cut in mm
V = 240/60; % cutting speed in m/sec

T = readtable('Book1.xlsx'); % table containing experimental data 
h = table2array(T(:,1));
Ft = table2array(T(:,2));
Ff = table2array(T(:,3));
hc = table2array(T(:,4));
figure(1);
plot(h,Ft,'b*'); hold on;
plot(h,Ff,'r*');


%% Iterating through feed values
% Initializing the arrays
phi_c = zeros(7,1);
beta_a = zeros(7,1);
F = zeros(7,1);
tau_s = zeros(7,1);
sigma_s = zeros(7,1);
Pt = zeros(7,1);
Ktc = zeros(7,1);
Kfc = zeros(7,1);
gamma_s = zeros(7,1); 
gamma_r = zeros(7,1); 
r_c = zeros(7,1); 

for i=1:7
    rc = h(i)/hc(i); % Chip thickness ratio
    r_c(i) = rc;
    phi_c(i) = atan(rc*cos(alpha_r)/(1-rc*sin(alpha_r))); % shear angle
    beta_a(i) = alpha_r + atan(Ff(i)/Ft(i)); % Friction angle
    F(i) = sqrt(Ft(i)^2 + Ff(i)^2); % Resultant Cutting Force
    Fs = F(i)*cos(phi_c(i) + beta_a(i) - alpha_r); % Shearing Force
    tau_s(i) = Fs/(b*h(i)/sin(phi_c(i))); % Shearing stress
    Fn = F(i)*sin(phi_c(i) + beta_a(i) - alpha_r); % Normal Force on shearing Plane
    sigma_s(i) = Fn/(b*h(i)/sin(phi_c(i))); % Normal stress on shearing Plane
    Vs = V*cos(alpha_r)/cos(phi_c(i)-alpha_r); % Shearing Velocity
    Ps = Fs*Vs; % Shearing Power
    Pu = F(i)*sin(beta_a(i))*(rc*V); % Frictional Power
    Pt(i) = Ps + Pu; % Total Cutting Power
    Ktc(i) = Ft(i)/(b*h(i)); % Cutting Force Coefficients
    Kfc(i) = Ff(i)/(b*h(i));
    gamma_s(i) = cos(alpha_r)/(sin(phi_c(i))*cos(phi_c(i)-alpha_r)); % Shear strain
    Lc = hc(i)/cos(phi_c(i)-alpha_r); % shaer plane length
    del_d = 0.15*Lc; % shear zone thickness
    gamma_r(i) = V*cos(alpha_r)/(del_d*cos(phi_c(i)-alpha_r)); % shear strain rate
end

Tfin = table(h,Ft,Ff,hc,phi_c, beta_a, F, tau_s, sigma_s, Pt, Ktc, Kfc, gamma_s, gamma_r); % Final Table as asked in Part 1

%% Curve Fitting
% Linear Model
mdFt = polyfit(h, Ft, 1); % linear regression model fit to Ft
mdFf = polyfit(h, Ff, 1); % linear regression model fit to Ff
xrange = linspace(0.02,0.08);
Ft_lin = polyval(mdFt, xrange);
Ff_lin = polyval(mdFf, xrange);
plot(xrange, Ft_lin, 'b');
plot(xrange, Ff_lin, 'r');
Ktc_lin = mdFt(1)/b;
Kte_lin = mdFt(2)/b;
Kfc_lin = mdFf(1)/b;
Kfe_lin = mdFf(2)/b;

% Non-linear Model
Ktc_fit = fit(h,Ktc,'exp1');
Kfc_fit = fit(h,Kfc,'exp1');
plot(xrange,Ktc_fit.a*exp(Ktc_fit.b*xrange)*b.*xrange,'b--')
plot(xrange,Kfc_fit.a*exp(Kfc_fit.b*xrange)*b.*xrange,'r--')

grid on;
lgd = legend('Experimental Values of F_t','Experimental Values of F_f', ... 
    'Linear Model: F_t','Linear Model: F_f','Non-Linear Model: F_t','Non-Linear Model: F_f');
lgd.NumColumns = 2;
lgd.Orientation = 'horizontal';
lgd.Location = 'northwest';
xlabel('Feed Rate in mm/rev');
ylabel('Forces in N');
title('Plot showing experimental values of tangential and feed forces in comparison to linear and non-linear models');

% Comparing the force coefficients
figure(2)
subplot(1,2,1)
plot(h,Ktc,'r*'); hold on;
plot(xrange, Ktc_fit.a*exp(Ktc_fit.b*xrange), 'r')
plot(xrange,Ktc_lin, 'r.')
grid on;
xlabel('Feed Rate in mm/rev');
ylabel('K_{tc} in N/mm^2');
legend('From empirical formula','From non-linear force model','From linear force model')
title('Plot showing comparison between force coefficient: K_{tc}');

subplot(1,2,2)
plot(h,Kfc,'b*'); hold on;
plot(xrange, Kfc_fit.a*exp(Kfc_fit.b*xrange), 'b')
plot(xrange,Kfc_lin, 'b.')
grid on;
xlabel('Feed Rate in mm/rev');
ylabel('K_{fc} in N/mm^2');
legend('From empirical formula','From non-linear force model','From linear force model')
title('Plot showing comparison between force coefficient: K_{fc}');

%% Shear angle comparison
phi_c_merch = pi/4 + alpha_r/2 - beta_a/2;
figure(3);
plot(h,phi_c_merch,'rd-'); hold on;
plot(h,phi_c,'bo-'); hold on;
grid on;
xlabel('Feed Rate in mm/rev');
ylabel('Shear angle \phi_{c}');
legend('Theoretical prediction of shear angle','Shear angle from orthogonal cutting tests','Location','northwest')
title('Shear angle comparison')

% figure(4)
% subplot(1,3,1)
% plot(h,tau_s,'b*'); hold on;
% grid on;
% xlabel('Feed Rate in mm/rev');
% ylabel('Shear stress \tau_{s}');
% mdtau_s = polyfit(h,tau_s,1);
% tau_s_lin = polyval(mdtau_s, xrange);
% plot(xrange,tau_s_lin,'b--');
% 
% subplot(1,3,2)
% plot(h,beta_a,'b'); hold on;
% grid on;
% xlabel('Feed Rate in mm/rev');
% ylabel('Friction angle \beta_{a}');
% mdbeta_a = polyfit(h,beta_a,1);
% beta_a_lin = polyval(mdbeta_a, xrange);
% plot(xrange,beta_a_lin,'b--');
% 
% subplot(1,3,3)
% plot(h,r_c,'b'); hold on;
% grid on;
% xlabel('Feed Rate in mm/rev');
% ylabel('Chip thickness ratio \r_{c}');
% mdr_c = polyfit(h,r_c,1);
% r_c_lin = polyval(mdr_c, xrange);
% plot(xrange,r_c_lin,'b--');
