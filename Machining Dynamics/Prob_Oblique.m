clc;
clear;

alpha_n = 20*(pi/180);
beta_a = 32.3525*(pi/180);

v = 0.1;
for ang = 0:5:30
    i = ang*(pi/180);
    eta_temp = i;
    for k = 0:10
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
        %plot(k,Pt_p,'bo'); hold on;
    end
    subplot(2,2,1)
    plot(ang,theta_n_val*(180/pi),'b.'); hold on;
    subplot(2,2,2)
    plot(ang,theta_i_val*(180/pi),'b.'); hold on;
    subplot(2,2,3)
    plot(ang,eta_temp*(180/pi),'b.'); hold on;
    subplot(2,2,4)
    plot(ang,phi_i_val*(180/pi),'b.'); hold on;    
end
 