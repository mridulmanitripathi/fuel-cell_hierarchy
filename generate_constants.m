function p = generate_constants()

p = struct; % paramter struct

%Constants
p.Ha     = 1e-5;        %Length of primary pores, cm     (50-100 nm)
p.L      = 1e-3;        %Length of secondary pores (= CCL thickness), cm
p.Ra     = 50*1e-7;     %Agglomerate radius Ra, cm
p.sigmaP = 0.0221;      %CCL proton conductivity Ohm cm-1
p.b      = 0.0302;      %ORR Tafel slope (V)
p.Cdl    = 20;          %F cm-3 from Kulikovsky12
p.F      = 96485.33212; %Faraday constant, C mol−1
p.Da     = 4*1e-5;      %Oxygen diffusion coefficient in agglomerate (in water) cm2 s−1
p.Dox    = 1e-4;        %Oxygen diffusion coefficient through the CCL depth, cm2 s−1 (upper estimate)
p.i_star = 1e-2;        %Exchange current density i∗, A cm−3

% p.p = 101300; %[N/m^2=Pa]
% p.T = 25+273.15; %[K]

p.c_ref =1e-6;         %reference molar concentration (mol cm-3) 10^-6 (i think from kulikovsky 15)

p.t_star = 4*p.F*p.c_ref/p.i_star;
p.C1 = p.Ha^2/p.L^2;
p.C2 = -3*p.Ha/(p.c_ref*p.Ra);
p.C3 = p.sigmaP*p.b/(p.Cdl*p.L^2)*p.t_star;
p.C4 = 12*(p.F*p.Ha)/(p.Ra*p.Cdl*p.b);
p.D_star = p.Ha^2*p.i_star/(4*p.F*p.c_ref^2);
p.Da_tilde = p.Da/p.D_star;
p.Dox_tilde = p.Dox/p.D_star;

end
