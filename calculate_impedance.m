%% Symbolic variables

syms ca(h) cox(x) eta(x) omega Da Dox b eta0 C C1 C2 C3 C4 t_star ca_cox J

% J := Dca(0)
% ca_cox := cox(x)

Dca  = diff(ca,h);
Dcox = diff(cox,x);
Deta = diff(eta,x);

%% Problem
eqn1 = (1i*omega*ca - Da*diff(ca,h,2)       == 0);
bv1  = [ca(0) == ca_cox
        Dca(1) == 0];

eqn2 = (1i*omega*cox - C1*Dox*diff(cox,x,2) == C2*Da*J);
bv2  = [Dcox(0) == 0
        cox(1) == 1];

eqn3 = (1i*omega*eta - C3*diff(eta,x,2)     == C4*Da*J);
bv3  = [eta(0) == eta0/b
        Deta(1) == 0];

%% Try to solve independently
ca_sol(h)  = dsolve(eqn1,bv1);
cox_sol(x) = dsolve(eqn2,bv2);
eta_sol(x) = dsolve(eqn3,bv3);

%% use definitions of the parameters and substitute them into the solution

% use J == Dca(0)
Dca_sol = diff(ca_sol, h);
J_sol = Dca_sol(0); % J depends only on ca_cox

% use ca_cox == cox_sol(0)  %TODO: IN ORDER TO OBTAIN A FUNCTION FOR ALL PRIMARY PORES, THIS SHOULD BE DONE FOR THE CORRESPONDING VALUES OF x
cox_sol_0_subs = subs(cox_sol(0), J, J_sol);
ca_cox_sol = solve(ca_cox == cox_sol_0_subs, ca_cox);

% substitute into ca_sol, cox_sol, eta_sol:
subs_fcn = @(expr) subs(subs(expr, J, J_sol), ca_cox, ca_cox_sol);

%% final system solution
ca_fin  = subs_fcn(ca_sol);
cox_fin = subs_fcn(cox_sol);
eta_fin = subs_fcn(eta_sol);

%% we need to get eta(0) and Deta(0) for computing the impedance
Deta_fin   = diff(eta_fin, x);
eta_fin_0  = eta_fin(0);
Deta_fin_0 = Deta_fin(0);

%% check if the solution is correct
% preparation:
Dca_fin = diff(ca_fin, h);
J_fin = Dca_fin(0);

% now check it:
% subs(eqn1, ca, ca_fin)
% subs(subs(eqn2, cox, cox_fin), J, J_fin)
% subs(subs(eqn3, eta, eta_fin), J, J_fin)
% subs(bv1,ca,ca_fin)
% subs(bv2,cox,cox_fin)
% subs(bv3,eta,eta_fin)

%% calculate the impedance
Z_tilde = -C*eta_fin_0/Deta_fin_0;

disp(Z_tilde);

%% get constants
p = generate_constants();

%% set parameters
p.eta0  = 0.4; % not sure
p.C     = 1; % 1e-5 (eqn. (14) in Kulikovsky17)

%% evaluation
close all;

%% 2) plot oxygen concentration in secondary pores
p.omega = 2*pi; % for example
p.omega_tilde = p.omega/p.t_star;
p.x_tilde = linspace(0,1,1000);

cox_tilde_subs = subs(cox_fin, [  Da         Dox         C1   C2], ...
                               [p.Da_tilde p.Dox_tilde p.C1 p.C2]);
cox_tilde_fcn = matlabFunction(cox_tilde_subs);
cox_tilde_fin = cox_tilde_fcn(p.x_tilde, p.omega_tilde);

figure(2);
plot(p.x_tilde, cox_tilde_fin);
title('Oxygen concentration in secondary pores');
xlabel('x tilde');
ylabel('cox tilde');
grid on;

%% 3) plot overpotential (in secondary pores)
p.omega = 2*pi; % for example
p.omega_tilde = p.omega/p.t_star;
p.x_tilde = linspace(0,1,1000);

eta_tilde_subs = subs(eta_fin, [  Da         Dox         b   C1   C2   C3   C4   t_star   eta0], ...
                               [p.Da_tilde p.Dox_tilde p.b p.C1 p.C2 p.C3 p.C4 p.t_star p.eta0]);
eta_tilde_fcn = matlabFunction(eta_tilde_subs);
eta_tilde_fin = eta_tilde_fcn(p.x_tilde, p.omega_tilde);

figure(3);
plot(p.x_tilde, eta_tilde_fin);
title('Overpotential (in secondary pores)');
xlabel('x tilde');
ylabel('eta tilde');
grid on;

%% 4) plot impedance
p.f = 10.^(linspace(0,3,1e6)); % (1e+0:1e+0:1e+6); 
p.omega = 2*pi*p.f;
p.omega_tilde = p.omega/p.t_star;

Z_tilde_subs = subs(Z_tilde, [  Da         Dox         b   C1   C2   C3   C4   t_star   eta0   C], ...
                             [p.Da_tilde p.Dox_tilde p.b p.C1 p.C2 p.C3 p.C4 p.t_star p.eta0 p.C]);
Z_tilde_fcn = matlabFunction(Z_tilde_subs);
Z_tilde_fin = Z_tilde_fcn(p.omega_tilde);

figure(4);
plot(real(Z_tilde_fin), -imag(Z_tilde_fin));
% plot3(real(Z_tilde_fin), -imag(Z_tilde_fin), p.omega);
title(sprintf('Impedance Z tilde for frequency f in [%d, %d] ', p.f(1), p.f(end)));
xlabel('Re(Z)');
ylabel('-Im(Z)');
zlabel('omega');
grid on;
