function dydt = MacArthurRosenzweig_bdLogistic_ExplicitTraits(~, y, b, d, qb, qd, xa, amax, theta_a, tau_a, e, xh, hmin, theta_h, tau_h, m)
dydt = zeros(size(y));

% variables

R = y(1);
C = y(2);
a = amax*exp(-((xa - theta_a)^2)/(2*tau_a^2));
h = hmin*exp(((xh - theta_h)^2)/(2*tau_h^2));

dydt(1) = R*((b - qb*R) - (d + qd*R)) - a*R*C/(1 + a*h*R);

dydt(2) = C*(e*a*R/(1 + a*h*R) - m);
