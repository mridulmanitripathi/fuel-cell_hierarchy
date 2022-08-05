clear;
close;

%% User input
approach = 'bvp4c';
ca_cox = 1;

%% Solve the system
switch approach
    %% Initial value problem (IVP) approach
    case 'ode45'
        [h,y] = ode45(@fuel, [1 0], [ca_cox; 0]);
        % Note: 

    %% Boundary value problem (BVP) approach
    case 'bvp4c'
        bv_ox_a = @(y0, y1, ca_cox) [y1(2); y0(1)-ca_cox]; % boundary value
        solinit = bvpinit(linspace(0,1,100), [ca_cox, 0]); % initial guess
        sol = bvp4c(@fuel, @(y0,y1) bv_ox_a(y0,y1,ca_cox), solinit);
        h = sol.x';
        y = sol.y';
end

%% Plot
plot(h,real(y(:,1)),'-',h,real(y(:,2)),'-');
title('Solution of fuel cell');
xlabel('length h');
ylabel('Solution c_a');
legend('c_a','d c_a');
