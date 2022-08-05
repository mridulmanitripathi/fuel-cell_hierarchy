%{ 
  Iterative Solver for system of cuopled ODEs using ode45 for secondary 
  pores and analytical solution for primary pores
%}
%% Define constants

p = generate_constants;

flux_pp = -0.05;
omega = 1000;
omega_tilde = omega*p.i_star;
Da = p.Da_tilde;
C1 = p.C1;
C2 = p.C2;
C3 = p.C3;
C4 = p.C4;
ca_cox = 1;

%Flux calculation from analytical solution
Flux_PP = - (ca_cox*(Da*omega_tilde*1i)^(1/2))/Da + (2*ca_cox*(Da*omega_tilde*1i)^(1/2))/(Da*(exp((2*(Da*omega_tilde*1i)^(1/2))/Da) + 1));
C_0 = (C2)*Da*Flux_PP;
%% Equation

syms cox(x)

[V] = odeToVectorField(diff(cox,2) == (-C_0 + 1i*omega_tilde*cox));
M = matlabFunction(V,'vars', {'x','Y'});



%% Solve in the interval [x_i , 1]

x_i = 0;
x_f = 1;

% Initial Conditions

cox_i = 1;           % local oxygen concentration
d2cox_i = 0.1;       % guesstimate local second derivative
cox_f = 0.1;   

cox_sol_i = ode45(M,[x_i x_f],[cox_i d2cox_i]);

range_1 = linspace(0,1,30);
cox_i = deval(cox_sol_i,range_1,1);

y_1 =(real(cox_i));
plot(range_1, y_1)

%% discretization
n = 50; % number of steps
step = 1/n;  %increment
%% Initialization
ca_cox = 1;
eta_x = 0.2;
concentration = zeros(1,n);
overpotential = zeros(1,n);
concentration(1)= ca_cox;
overpotential(1) = eta_x;
index_1 = 1;
x_f =1.1;
vec1= 0.02:step:0.98;
%% Solve iteratively 
for it = vec1
    %compute flux at the boundary with primary pores for position x_i
    Flux_PP = - (ca_cox*(Da*omega_tilde*1i)^(1/2))/Da + (2*ca_cox*(Da*omega_tilde*1i)^(1/2))/(Da*(exp((2*(Da*omega_tilde*1i)^(1/2))/Da) + 1));
    C_0 = (C2)*Da*Flux_PP;
    C_01 = (C4)*Flux_PP;
    
    % interval for solution
    x_i = it;
    
    %% Equation for oxygen concentration
    syms cox(x)
    [V] = odeToVectorField(diff(cox,2) == (-C_0 + 1i*omega_tilde*cox));
    M = matlabFunction(V,'vars', {'x','Y'});
    % Initial conditions
    cox_i = ca_cox;           % local oxygen concentration
    d2cox_i = 0.1;            % guesstimate local second derivative
    
    %Solve the ODE
    cox_sol_i = ode45(M,[x_i x_f],[cox_i d2cox_i]);
    %evaluate at the next point
    eval_point = ((1-it)*0.07)+it;
    ca_cox = real(deval(cox_sol_i,eval_point,1));
    index_1 = index_1 + 1;
    concentration(index_1)=ca_cox;
    
    %% Equation for overpotential
    eta_i = eta_x;
    syms eta(x)
    [V] = odeToVectorField(diff(eta,2) == (C_01+(1i*omega_tilde*eta)));
    N = matlabFunction(V,'vars', {'x','Y'});
    %Solve the ODE
    eta_sol_i = ode45(N,[x_i x_f],[eta_i  (-d2cox_i/1.5)]);
    % Evaluation
    eta_x = real(deval(eta_sol_i,eval_point,1));
    % Store value
    overpotential(index_1)= eta_x;
    
    
end
%% Plot Results
x_i = linspace(0,1,n);
plot(x_i, concentration,x_i, overpotential)
title('Oxygen and Overpotential in CCL')
legend('C_{ox}','\eta(x)', 'Location', 'Northwest')
