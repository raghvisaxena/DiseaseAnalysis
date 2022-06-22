beta = 0.5; % Transmission Rate
gamma  = 0.14; % Recovery Rate
lambda = 0.5; % Rate at which exposed individuals become Infectious
N = 309;                          % Total population 
I0 = 29;                              % Initial number of infected individuals
E0 = 0;                               % Initial number of exposed individuals
S0 = N;                               % Initial number of susceptible individuals
T = 70;                              % Period of 70 days
R0=beta/gamma;                        % Reproduction number
% Solve the ODE using ode45.
y0 = [S0; E0; I0; R0];                % Initial conditions
tspan = [0 T];                        % Interval of Integration
[t,y] = ode45(@(t,y) seir_model(t,y,beta,lambda,gamma,N),tspan,y0);

% Plots that show the results
plot(t,y,'LineWidth',2);
xlabel('Days'); 
ylabel('Number of individuals');
legend('Suseptible','Exposed','Infected','Recovered');
title('Hyderabad-2018');

function dydt = seir_model(t,y,beta,lambda,gamma,N)
  S = y(1);
  E = y(2);
  I = y(3);
  % Equations of the model described above
  dS = (-beta*S*I/N) ;
  dE = (beta*I*S/N - lambda*E) ;
  dI = (lambda*E - gamma*I);
  dR = (gamma*I);
  dydt = [dS;dE;dI;dR];
end