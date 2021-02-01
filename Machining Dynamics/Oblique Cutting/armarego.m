clc;
clear;

alpha_n = 20*(pi/180);
beta_a = 32.3525*(pi/180);
h = 0.5; % given
hc = 1.1; % assumed
rc = h/hc;
v = 0.1;
for ang = 1:30
    i = ang*(pi/180);
    eta_temp_a = i;
%     syms eta phi_n beta_n
%     eqn1 = tan(phi_n + beta_n) == (cos(alpha_n)*tan(i))/(tan(eta) - sin(alpha_n)*tan(i));
%     eqn2 = tan(beta_n) == tan(beta_a)*cos(eta);
%     eqn3 = tan(phi_n) == (rc*(cos(eta)/cos(i))*cos(alpha_n))/(1-rc*(cos(eta)/cos(i))*sin(alpha_n));
    
    %SSA = vpasolve([eqn1, eqn2, eqn3], [eta phi_n beta_n], [0; 40*pi/180 ;0]);
    for k=1:10
        eta_a = eta_temp_a;
        beta_n_a = atan(tan(beta_a)*cos(eta_a));
        phi_n_a = atan((rc*(cos(eta_a)/cos(i))*cos(alpha_n))/(1-rc*(cos(eta_a)/cos(i))*sin(alpha_n)));
        if k ==1
            eta_s = eta_temp_a;
            phi_n_s = phi_n_a;
            beta_n_s = beta_n_a;
        end
        eta_temp_a = atan((cos(alpha_n)*tan(i))/(tan(phi_n_a + beta_n_a)) + (sin(alpha_n)*tan(i)));
        %plot(k,eta_temp,'b.'); hold on;
    end
    
    subplot(2,2,1)
    plot(ang,eta_temp_a*(180/pi),'b.'); hold on;
    plot(ang,eta_s*(180/pi),'r.'); hold on; % armarego simplified
    subplot(2,2,2)
    plot(ang,phi_n_a*(180/pi),'b.'); hold on;
    plot(ang,phi_n_s*(180/pi),'r.'); hold on;
    subplot(2,2,3)
    plot(ang,(beta_n_a-alpha_n)*(180/pi),'b.'); hold on;
    plot(ang,(beta_n_s-alpha_n)*(180/pi),'r.'); hold on;
end