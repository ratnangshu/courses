clc;
clear;

%% Experimental Data
T = readtable('Book1.xlsx'); % table containing experimental data 
Feed_arr = table2array(T(:,1));
Ft_arr = table2array(T(:,2));
Ff_arr = table2array(T(:,3));
Chipt_arr = table2array(T(:,4));
figure(1);
plot(Feed_arr,Ft_arr,'b*'); hold on;
plot(Feed_arr,Ff_arr,'r*');

%% Evaluate the cutting coefficients by a linear regression of the measured forces
alpha_r = 0; % rake angle
b = 5; % width of cut in mm
V = 240; % cutting speed in m/min
rho = 7800; % Density of the steel that was cut in kg/m^3

mdFt = polyfit(Feed_arr, Ft_arr, 1); % linear regression model fit to Ft
mdFf = polyfit(Feed_arr, Ff_arr, 1); % linear regression model fit to Ff

xrange = linspace(0,0.08);
Ft = polyval(mdFt, xrange);
Ff = polyval(mdFf, xrange);
plot(xrange, Ft, 'b');
plot(xrange, Ff, 'r');

Ktc = mdFt(1)/b;
Kte = mdFt(2)/b;
Kfc = mdFf(1)/b;
Kfe = mdFf(2)/b;

%% Evaluate the shear angle, shear stress, and avg friction coeff. for each test
phi_c = zeros(7,1);
tau_s = zeros(7,1);
mu_a = zeros(7,1);
gamma_s = zeros(7,1); % shear strain
gamma_r = zeros(7,1); % strain rate
del_d = 0.2; % shear zone thickness
for i=1:7
    h = Feed_arr(i);
    hc = Chipt_arr(i);
    rc = h/hc;
    phi_c(i) = atan(rc*cos(alpha_r)/(1-rc*sin(alpha_r)));
    tau_s(i) = (Ft_arr(i)*cos(phi_c(i)) - Ff_arr(i)*sin(phi_c(i)))/(b*h/sin(phi_c(i)));
    mu_a(i) = (Ft_arr(i)*sin(alpha_r) + Ff_arr(i)*cos(alpha_r))/(Ft_arr(i)*cos(alpha_r) - Ff_arr(i)*sin(alpha_r));
    
    gamma_s(i) = cot(phi_c(i)) + tan(phi_c(i) - alpha_r);
    gamma_r(i) = V*cos(alpha_r)/(del_d*cos(phi_c(i)-alpha_r));
end

%% Express them as an empirical function of uncut chip thickness (h)
mdphi_c = polyfit(Feed_arr, phi_c, 1); % linear regression model fit to phi_c
mdtau_s = polyfit(Feed_arr, tau_s, 1); % linear regression model fit to tau_s
mdmu_a = polyfit(Feed_arr, mu_a, 1); % linear regression model fit to mu_a
phi_c_fit = polyval(mdphi_c, xrange);
tau_s_fit = polyval(mdtau_s, xrange);
mu_a_fit = polyval(mdmu_a, xrange);

figure(2);
subplot(1,3,1)
plot(Feed_arr,phi_c,'rs'); hold on;
plot(xrange, phi_c_fit, 'r');
axis square
subplot(1,3,2)
plot(Feed_arr,tau_s,'gs'); hold on;
plot(xrange, tau_s_fit, 'g');
axis square
subplot(1,3,3)
plot(Feed_arr,mu_a,'bs'); hold on;
plot(xrange, mu_a_fit, 'r');
axis square

%% Predict the cutting force coefficients using empirically expressed formula

Ktc_emp = zeros(7,1);
Kfc_emp = zeros(7,1);
beta_a = zeros(7,1);
for i=1:7
    beta_a(i) = atan(mu_a(i));
    Ktc_emp(i) = tau_s(i)*(cos(beta_a(i)-alpha_r))/(sin(phi_c(i))*cos(phi_c(i)+beta_a(i)-alpha_r));
    Kfc_emp(i) = tau_s(i)*(sin(beta_a(i)-alpha_r))/(sin(phi_c(i))*cos(phi_c(i)+beta_a(i)-alpha_r));
end

%% Predicting shear from merchant
phi_c_merch = pi/4 + alpha_r/2 - beta_a/2;
figure(2);
subplot(1,3,1)
plot(Feed_arr,phi_c_merch,'rd'); hold on;


