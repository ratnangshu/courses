function [Ft, Ff, Fr] = calcForce(i, theta_i, theta_n, phi_i, phi_n, beta_n, alpha_n, eta, me)
ts = 350;
b = 6.25;
h = 0.5;
if me == 1 % me = 1 follow set of formulas for minimum energy method
    Ktc = ts*(cos(theta_n) + tan(theta_i)*tan(i))/((cos(theta_n+phi_n)*cos(phi_i) + tan(theta_i)*sin(phi_i))*sin(phi_n));
    Kfc = ts*(sin(theta_n))/((cos(theta_n+phi_n)*cos(phi_i) + tan(theta_i)*sin(phi_i))*sin(phi_n));
    Krc = ts*(tan(theta_i) - cos(theta_n)*tan(i))/((cos(theta_n+phi_n)*cos(phi_i) + tan(theta_i)*sin(phi_i))*-sin(phi_n));
else % me = 0 follow set of formulas for armarego method
    Ktc = (ts/sin(phi_n))*((cos(beta_n-alpha_n)+tan(i)*tan(eta)*tan(beta_n))/(sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)^2)));
    Kfc = (ts/(sin(phi_n)*cos(i)))*((sin(beta_n-alpha_n))/(sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)^2)));
    Krc = (ts/sin(phi_n))*((cos(beta_n-alpha_n)*cos(beta_n-alpha_n)*tan(i)-tan(eta)*tan(beta_n))/(sqrt(cos(phi_n+beta_n-alpha_n)^2 + tan(eta)^2*sin(beta_n)^2)));
end
Ft = Ktc*b*h;
Ff = Kfc*b*h;
Fr = Krc*b*h;
end

