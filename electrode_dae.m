function f = electrode_dae(X)

F = 96487;         % C/mol
R = 8.314;         % J/mol K
T = 298.15;        % K
phi_eq1 = 0.42;    % V
phi_eq2 = 0.303;   % V
rho = 3.4;         % g/cm^3
W = 92.7;          % g/mol
V = 1*10^(-5);  % cm
i_app = 1*10^(-5);  % A/cm^2
i01 = 1*10^(-4);    % A/cm^2
i02 = 1*10^(-8);    % A/cm^2

y1 = X(1);
y2 = X(2);



j1 = i01*(2*(1-y1)*exp((0.5*F/(R*T))*(y2-phi_eq1)) - 2*y1*exp((-0.5*F/(R*T))*(y2-phi_eq1)));
j2 = i02*(exp((F/(R*T))*(y2-phi_eq2)) - exp((-F/(R*T))*(y2-phi_eq2)));

f = [
       (W/(F*rho*V)) * j1;
       j1 + j2 - i_app;
     ];
 
end