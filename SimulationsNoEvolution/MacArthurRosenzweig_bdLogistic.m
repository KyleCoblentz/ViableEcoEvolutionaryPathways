function dydt = MacArthurRosenzweig_bdLogistic(~, y, b, d, qb, qd, a, e, h, m)
dydt = zeros(size(y));

% variables

R = y(1);
C = y(2);

dydt(1) = R*((b - qb*R) - (d + qd*R)) - a*R*C/(1 + a*h*R);

dydt(2) = C*(e*a*R/(1 + a*h*R) - m);


