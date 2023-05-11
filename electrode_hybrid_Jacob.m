function [A_x] = electrode_hybrid_Jacob( X_k ) 
 
electrode_global ; 
global wk vk w1 w2 w3


y1 = X_k(1);
y2 = X_k(2);
w1 = X_k(3);
w2 = X_k(4);
w3 = X_k(5);


dj1_by_dy1 = -2*i01*( exp(((0.5*F)/(R*T))*(y2 - phi_eq1)) + exp(((-0.5*F)/(R*T))*(y2 - phi_eq1)) );


dj1_by_dy2 = 2.*i01.*( ((1-y1).*((0.5.*F)/(R.*T)).*exp(((0.5.*F)/(R.*T)).*(y2 - phi_eq1))) + (y1.*((0.5.*F)/(R.*T)).*exp(((-0.5.*F)/(R.*T)).*(y2 - phi_eq1))) );
dj2_by_dy1 = 0;
dj2_by_dy2 = i02.*( ((F/(R.*T)).*exp((F/(R*T)).*(y2 - phi_eq2))) + ((F/(R.*T)).*exp(((-F)/(R.*T)).*(y2 - phi_eq2))) );

A1 = (W.*dj1_by_dy1)/(rho*V*F);
B1 = (W.*dj1_by_dy2)/(rho*V*F);
C1 = dj1_by_dy1 + dj2_by_dy1;
D1 = dj1_by_dy2 + dj2_by_dy2;


% A_x = [A1 B1 0 0 0;-D1^-1*C1*A1 -D1^-1*C1*B1 0 0 0; 0 0 w1 0 0; 0 0 0 w2 0; 0 0 0 0 w3];

A_x = [A1 B1 0 0 0;-D1^-1*C1*A1 -D1^-1*C1*B1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];

end