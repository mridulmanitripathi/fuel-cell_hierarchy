%{ 
  Impedance spectrum and Iterative Solver for system of cuopled ODEs using ode45 for secondary 
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





%% discretization
n = 50; % number of steps
step = 1/n;  %increment
n1 = 100;
impedance = zeros(n1,1);
imp_index =1;
for omega = 0.01:10:1000
    
    %% Initialization
    ca_cox = 1;
    eta_x = 0.1;
    concentration = zeros(1,n);
    overpotential = zeros(1,n);
    concentration(1)= ca_cox;
    overpotential(1) = eta_x;
    index_1 = 1;
    x_f =1.1;
    vec1= 0.02:step:0.98;
    omega_tilde = omega*p.i_star;
    %% Solve Equations
    
    %compute flux at the boundary with primary pores for position x_i
    Flux_PP = - (ca_cox*(Da*omega_tilde*1i)^(1/2))/Da + (2*ca_cox*(Da*omega_tilde*1i)^(1/2))/(Da*(exp((2*(Da*omega_tilde*1i)^(1/2))/Da) + 1));
    C_0 = (C2)*Da*Flux_PP;
    C_01 = (C4)*Flux_PP;

    % interval for solution
    x_i = 0;

    %% Equation for oxygen concentration
    syms cox(x)
    [V] = odeToVectorField(diff(cox,2) == (-C_0 + 1i*omega_tilde*cox)/C1);
    M = matlabFunction(V,'vars', {'x','Y'});
    % Initial conditions
    cox_i = ca_cox;           % local oxygen concentration
    dcox_i = 0.1;       % guesstimate local first derivative

    %Solve the ODE
    cox_sol_i = ode45(M,[x_i x_f],[cox_i dcox_i]);
    %evaluate at the next point
    eval_point = 0.3;
    ca_cox = real(deval(cox_sol_i,eval_point,1));
    index_1 = index_1 + 1;
    concentration(index_1)=ca_cox;

    %% Equation for overpotential
    eta_i = eta_x;
    syms eta(x)
    [V] = odeToVectorField(diff(eta,2) == (C_01+(1i*omega_tilde*eta))/C3);
    N = matlabFunction(V,'vars', {'x','Y'});
    %Solve the ODE
    eta_sol_i = ode45(N,[x_i x_f],[eta_i  (-dcox_i/1.5)]);
    % Evaluation
    eta_x = real(deval(eta_sol_i,eval_point,1));
    % Store value
    overpotential(index_1)= eta_x;
        



    eta_0 = deval(eta_sol_i,0.01,1);
    eta_1 = deval(eta_sol_i,0.02,1);
    deta = (eta_1-eta_0)/(0.01);
    imp_i = eta_0/(C3*deta_0);
    impedance(imp_index) = -imp_i;
    imp_index = imp_index +1;
    
end
%% Plot Results

plot(real(impedance), -imag(impedance))
title('Impedance')

