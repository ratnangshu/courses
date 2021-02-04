clc;
clear;

%% Values according to Lin & Oxley (1972)
alpha_n = 20*(pi/180); % normal rake angle
beta_a = 32.3525*(pi/180); % friction angle 
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

%% Start of loop
for ang = 1:1:45
    i = ang*(pi/180); % inclination angle
    eta_temp = i; % chip flow angle for "minimum energy" calculation part
    eta_temp_a = i; % chip flow angle for "armarego" calculation part
    for k = 0:30 % the iteration loop typically converges by 30 iterations; can be replaced by a while loop
%     k = 0;
%     while(abs(Pt_p - Pt_p_old)>0.002)
%         k = k+1;
%         Pt_p_old = Pt_p;
%         if k>45
%             break
%         end

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
%     %% Plotting
%     subplot(2,3,1)
%     plot(ang,theta_n_val*(180/pi),'b.'); hold on;
%     plot(ang,(beta_n_a - alpha_n)*(180/pi),'g.'); hold on;
%     plot(ang,(beta_n_s - alpha_n)*(180/pi),'r.'); hold on;
%     ylim([-20 20])
%     subplot(2,3,2)
%     plot(ang,theta_i_val*(180/pi),'b.'); hold on;
%     ylim([0 50])
%     subplot(2,3,3)
%     plot(ang,phi_n_val*(180/pi),'b.'); hold on;
%     plot(ang,phi_n_a*(180/pi),'g.'); hold on;
%     plot(ang,phi_n_s*(180/pi),'r.'); hold on;
%     ylim([0 80])
%     subplot(2,3,4)
%     plot(ang,phi_i_val*(180/pi),'b.'); hold on;   
%     ylim([0 50])
%     subplot(2,3,5)
%     plot(ang,eta_temp*(180/pi),'b.'); hold on;
%     plot(ang,eta_temp_a*(180/pi),'g.'); hold on;
%     plot(ang,eta_s*(180/pi),'r.'); hold on; % armarego simplified
%     ylim([0 80])
end

ang_arr = 1:1:45;
%% Plotting
subplot(2,3,1)
plot(ang_arr,theta_n_val_arr*(180/pi),'b'); hold on;
plot(ang_arr,(beta_n_a_arr - alpha_n)*(180/pi),'m'); hold on;
plot(ang_arr,(beta_n_s_arr - alpha_n)*(180/pi),'r'); hold on;
xlim([0 45])
ylim([-20 20])
xlabel('Inclination angle i (degrees)') 
ylabel('\theta_n (degrees)') 
grid on

subplot(2,3,2)
plot(ang_arr,theta_i_val_arr*(180/pi),'b'); hold on;
xlim([0 45])
ylim([0 50])
xlabel('Inclination angle i (degrees)') 
ylabel('\theta_i (degrees)') 
grid on

subplot(2,3,4)
plot(ang_arr,phi_n_val_arr*(180/pi),'b'); hold on;
plot(ang_arr,phi_n_a_arr*(180/pi),'m'); hold on;
plot(ang_arr,phi_n_s_arr*(180/pi),'r'); hold on;
xlim([0 45])
ylim([0 80])
xlabel('Inclination angle i (degrees)') 
ylabel('\phi_n (degrees)') 
grid on

subplot(2,3,5)
plot(ang_arr,phi_i_val_arr*(180/pi),'b'); hold on;  
xlim([0 45])
ylim([0 50])
xlabel('Inclination angle i (degrees)') 
ylabel('\phi_i (degrees)') 
grid on

subplot(2,3,6)
plot(ang_arr,eta_temp_arr*(180/pi),'b'); hold on;
plot(ang_arr,eta_temp_a_arr*(180/pi),'m'); hold on;
plot(ang_arr,eta_s_arr*(180/pi),'r'); hold on; % armarego simplified
xlim([0 45])
ylim([0 80])
xlabel('Inclination angle i (degrees)') 
ylabel('\eta (degrees)') 
grid on

%% Plot Experimental Values
T = readtable('Book1.xlsx');
irange = 0:10:30;
subplot(2,3,6)% plot eta
etarange = [0 9.5 20 28.7];
plot(irange,etarange,'k*')
% plot phi_n
subplot(2,3,4)
phi_n_range = [26.3 27.2 27.4 27];
plot(irange,phi_n_range,'k*')

set(gcf, 'Position',  [180, 240, 1280, 360])

annotation('textbox',[0.7 0.575 .15 .4],'BackgroundColor','w');
annotation('textbox',[0.7 0.65 .3 .3],'String','\bullet Minimum energy principle','FitBoxToText','on','EdgeColor','none','color','b');
annotation('textbox',[0.7 0.55 .3 .3],'String','\bullet Armarego Model','FitBoxToText','on','EdgeColor','none','color','r');
annotation('textbox',[0.7 0.45 .3 .3],'String','\bullet Simplified Armarego Model','FitBoxToText','on','EdgeColor','none','color','m');
annotation('textbox',[0.7 0.35 .3 .3],'String','* Lin & Oxley','FitBoxToText','on','EdgeColor','none','color','k');