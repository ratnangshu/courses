clc;
clear;

%% Values according to Lin & Oxley (1972)
alpha_n = 20*(pi/180); % normal rake angle
beta_a = 32.57375*(pi/180); % friction angle, taken average
h = 0.5; % given depth of cut
hc = 1.1; % assumed cut chip thickness
rc = h/hc; % chip thickness ratio
v = 0.1; % interpolation ratio selected within the range 0 < v < 1

%% Initializing arrays
theta_n_val_arr = zeros(1,45);
beta_n_a_arr = zeros(1,45);
beta_n_s_arr = zeros(1,45);
theta_i_val_arr = zeros(1,45);
phi_n_val_arr = zeros(1,45);
phi_n_a_arr = zeros(1,45);
phi_n_s_arr = zeros(1,45);
phi_i_val_arr = zeros(1,45);
eta_temp_arr = zeros(1,45);
eta_temp_a_arr = zeros(1,45);
eta_s_arr = zeros(1,45);
Ft_arr = zeros(1,45);
Ff_arr = zeros(1,45);
Fr_arr = zeros(1,45);
Ft_a_arr = zeros(1,45);
Ff_a_arr = zeros(1,45);
Fr_a_arr = zeros(1,45);
Ft_s_arr = zeros(1,45);
Ff_s_arr = zeros(1,45);
Fr_s_arr = zeros(1,45);

%% Start of loop
for ang = 1:1:45
    i = ang*(pi/180); % inclination angle
    eta_temp = i; % chip flow angle for "minimum energy" calculation part
    eta_temp_a = i; % chip flow angle for "armarego" calculation part
    for k = 0:30 % the iteration loop typically converges by 30 iterations; can be replaced by a while loop
        
        %% Minimum Energy Principle
        % Force Relation
        eta_val = eta_temp;
        theta_i_val = asin(sin(beta_a)*sin(eta_val));
        theta_n_val_temp = atan(cos(eta_val)*tan(beta_a)) - alpha_n;
        theta_n_val = atan(sin(theta_n_val_temp)/cos(theta_n_val_temp));

        % Shear angle prediction
        syms phi_n phi_i
        Pt_p = (cos(theta_n_val) + tan(theta_i_val)*tan(i))/((cos(theta_n_val + phi_n)*cos(phi_i) + tan(theta_i_val)*sin(phi_i))*sin(phi_n));
        eqn4 = diff(Pt_p, phi_n) == 0;
        eqn5 = diff(Pt_p, phi_i) == 0;
        SSA = vpasolve([eqn4, eqn5], [phi_i, phi_n], [-pi/2 pi/2; -pi/2 pi/2]);
        %SSA = vpasolve([eqn4, eqn5], [phi_i, phi_n], [0; 40*pi/180]);
        phi_i_val = double(SSA.phi_i);
        phi_n_val = double(SSA.phi_n);

        % Velocity relation
        eta_new = atan(double(subs((tan(i)*cos(phi_n_val - alpha_n) - cos(alpha_n)*tan(phi_i_val))/(sin(phi_n_val)))));
        eta_temp = v*eta_temp + (1-v)*eta_new;
        
        % Power 
        Pt_p = (cos(theta_n_val) + tan(theta_i_val)*tan(i))/((cos(theta_n_val + phi_n_val)*cos(phi_i_val) + tan(theta_i_val)*sin(phi_i_val))*sin(phi_n_val));
        % plot(k,Pt_p,'bo'); hold on; % this could be plotted for separate
        % cases to see how power converges with the iterations
        
        %% Armarego 
        eta_a = eta_temp_a;
        beta_n_a = atan(tan(beta_a)*cos(eta_a));
        phi_n_a = atan((rc*(cos(eta_a)/cos(i))*cos(alpha_n))/(1-rc*(cos(eta_a)/cos(i))*sin(alpha_n)));
        if k ==1 % for simplified armarego case, which directly assumes eta = i;
            eta_s = eta_temp_a;
            phi_n_s = phi_n_a;
            beta_n_s = beta_n_a;
        end
        eta_temp_a = atan((cos(alpha_n)*tan(i))/(tan(phi_n_a + beta_n_a)) + (sin(alpha_n)*tan(i)));
    end
    
    %% Force Caluclations
    
    % Minimum Energy
    [Ft_arr(ang), Ff_arr(ang), Fr_arr(ang)] = calcForce(i, theta_i_val, theta_n_val, phi_i_val, phi_n_val, 0, alpha_n, eta_temp, 1);
    
    % Armarego
    [Ft_a_arr(ang), Ff_a_arr(ang), Fr_a_arr(ang)] = calcForce(i, 0, 0, 0, phi_n_a, beta_n_a, alpha_n, eta_temp_a, 0);
    
    % simplified Armarego
    [Ft_s_arr(ang), Ff_s_arr(ang), Fr_s_arr(ang)] = calcForce(i, 0, beta_n_s-alpha_n, 0, phi_n_s, beta_n_s, alpha_n, eta_s, 0);
    
    %% Storing values in arrays
    theta_n_val_arr(ang) = theta_n_val;
    beta_n_a_arr(ang) = beta_n_a;
    beta_n_s_arr(ang) = beta_n_s;
    theta_i_val_arr(ang) = theta_i_val;
    phi_n_val_arr(ang) = phi_n_val;
    phi_n_a_arr(ang) = phi_n_a;
    phi_n_s_arr(ang) = phi_n_s;
    phi_i_val_arr(ang) = phi_i_val;
    eta_temp_arr(ang) = eta_temp;
    eta_temp_a_arr(ang) = eta_temp_a;
    eta_s_arr(ang) = eta_s;
end

ang_arr = 1:1:45;
%% Plotting
subplot(3,3,1)
plot(ang_arr,theta_n_val_arr*(180/pi),'b'); hold on;
plot(ang_arr,(beta_n_a_arr - alpha_n)*(180/pi),'m'); hold on;
plot(ang_arr,(beta_n_s_arr - alpha_n)*(180/pi),'r'); hold on;
xlim([0 45])
ylim([-20 20])
xlabel('Inclination angle i (degrees)') 
ylabel('\theta_n (degrees)') 
grid on
grid minor

subplot(3,3,2)
plot(ang_arr,theta_i_val_arr*(180/pi),'b'); hold on;
xlim([0 45])
ylim([0 50])
xlabel('Inclination angle i (degrees)') 
ylabel('\theta_i (degrees)') 
grid on
grid minor

subplot(3,3,4)
plot(ang_arr,phi_n_val_arr*(180/pi),'b'); hold on;
plot(ang_arr,phi_n_a_arr*(180/pi),'m'); hold on;
plot(ang_arr,phi_n_s_arr*(180/pi),'r'); hold on;
xlim([0 45])
ylim([0 80])
xlabel('Inclination angle i (degrees)') 
ylabel('\phi_n (degrees)') 
grid on
grid minor

subplot(3,3,5)
plot(ang_arr,phi_i_val_arr*(180/pi),'b'); hold on;  
xlim([0 45])
ylim([0 50])
xlabel('Inclination angle i (degrees)') 
ylabel('\phi_i (degrees)') 
grid on
grid minor

subplot(3,3,6)
plot(ang_arr,eta_temp_arr*(180/pi),'b'); hold on;
plot(ang_arr,eta_temp_a_arr*(180/pi),'m'); hold on;
plot(ang_arr,eta_s_arr*(180/pi),'r'); hold on; % armarego simplified
xlim([0 45])
ylim([0 80])
xlabel('Inclination angle i (degrees)') 
ylabel('\eta (degrees)') 
grid on
grid minor

subplot(3,3,7)
plot(ang_arr, Ft_arr/1000, 'b'); hold on;
plot(ang_arr, Ft_a_arr/1000, 'm'); hold on;
plot(ang_arr, Ft_s_arr/1000, 'r'); hold on;
xlim([0 45])
ylim([0 5])
xlabel('Inclination angle i (degrees)') 
ylabel('F_{t} (kN)') 
grid on
grid minor

subplot(3,3,8)
plot(ang_arr, Ff_arr/1000, 'b'); hold on;
plot(ang_arr, Ff_a_arr/1000, 'm'); hold on;
plot(ang_arr, Ff_s_arr/1000, 'r'); hold on;
xlim([0 45])
ylim([0 1])
xlabel('Inclination angle i (degrees)') 
ylabel('F_{f} (kN)') 
grid on
grid minor

subplot(3,3,9)
plot(ang_arr, Fr_arr/1000, 'b'); hold on;
plot(ang_arr, Fr_a_arr/1000, 'm'); hold on;
plot(ang_arr, Fr_s_arr/1000, 'r'); hold on;
xlim([0 45])
ylim([0 5])
xlabel('Inclination angle i (degrees)') 
ylabel('F_{r} (kN)') 
grid on
grid minor

% Plotting Experimental Values
T = readtable('Book1.xlsx'); % table containing experimental data from lin and oxley 1972
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; 	[0.4940, 0.1840, 0.5560]; 	[0.4660, 0.6740, 0.1880]; 	[0.3010, 0.7450, 0.9330]];
irange = 0:10:30;
for sp_k = 3:4 % iterating for data of different cutting velocities U
    subplot(3,3,6) % plot eta data
    etarange = table2array(T(sp_k:7:sp_k+3*7,4));
    plot(irange,etarange,'*','color',colors(sp_k,:)); hold on;
    
    subplot(3,3,4) % plot phi_n data
    phi_n_range = table2array(T(sp_k:7:sp_k+3*7,6));
    plot(irange,phi_n_range,'*','color',colors(sp_k,:)); hold on;
    
    subplot(3,3,7) % plot F_t data
    F_t_range = table2array(T(sp_k:7:sp_k+3*7,10))/1000;
    plot(irange,F_t_range,'*','color',colors(sp_k,:)); hold on;
    
    subplot(3,3,8) % plot F_f data
    F_f_range = table2array(T(sp_k:7:sp_k+3*7,11))/1000;
    plot(irange,F_f_range,'*','color',colors(sp_k,:)); hold on;
    
    subplot(3,3,9) % plot F_r data
    F_r_range = table2array(T(sp_k:7:sp_k+3*7,12))/1000;
    plot(irange,F_r_range,'*','color',colors(sp_k,:)); hold on;
    
end

set(gcf, 'Position',  [180, 240, 1280, 360])

annotation('textbox',[0.65 0.7 .15 .25],'BackgroundColor','w');
annotation('textbox',[0.65 0.65 .3 .3],'String','\bullet Minimum Energy Principle','FitBoxToText','on','EdgeColor','none','color','b');
annotation('textbox',[0.65 0.6 .3 .3],'String','\bullet Armarego Model','FitBoxToText','on','EdgeColor','none','color','r');
annotation('textbox',[0.65 0.55 .3 .3],'String','\bullet Simplified Armarego Model','FitBoxToText','on','EdgeColor','none','color','m');
annotation('textbox',[0.65 0.5 .3 .3],'String','* Lin & Oxley (1972)','FitBoxToText','on','EdgeColor','none','color','k');

annotation('textbox',[0.825 0.8 .15 .15],'BackgroundColor','w');
annotation('textbox',[0.825 0.65 .3 .3],'String','* U = 60.96 m/min','FitBoxToText','on','EdgeColor','none','color',colors(2,:));
annotation('textbox',[0.825 0.6 .3 .3],'String','* U = 121.92 m/min','FitBoxToText','on','EdgeColor','none','color',colors(3,:));