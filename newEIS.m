clc
clear
lt=10^-3;
sigmap=0.02;
b=0.03;

n=1e5;
f=1:1e-2:1e4;
w=2*pi*f;
%Ha=50e-7;
ha = [50 60 70 80 90 100]'*(1e-7);  % 50-100nm
for tt=1:6
    Ha=ha(tt);
Ra = 1e-5;
C_ref = 1e-6;
i_star = 1e-2;
F = 96485.33;
Cdl = 20;

t_star = 4*F*C_ref/i_star;
C1 = (Ha^2)/(lt^2);
C2 = -3*(Ha)/Ra;
C3 = sigmap*b*t_star/(Cdl*lt^2);
D_star = i_star*(Ha^2)/(4*F*C_ref^2);
Da = (4e-5)/D_star;
Dox = (1e-4)/D_star;
eta0 = -0.42; %very less value gives more visual result
C4 = 12*F*Ha/(Cdl*b*Ra);
C=1;


for m=1:n
    z(m) = -(C*(eta0*w(m)*1i + eta0*w(m)*exp((2*(C3*w(m)*1i)^(1/2))/C3)*1i)*1i)/(b*w(m)*((((eta0*w(m)*(C3*w(m)*1i)^(1/2)*2i)/C3 + (C4*Da*b*((2*(Da*w(m)*1i)^(1/2))/(Da*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))*(1 + (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*2i)/(w(m)*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))) - (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*1i)/w(m))) - (4*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))*(1 + (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*2i)/(w(m)*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))) - (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*1i)/w(m))))*(C3*w(m)*1i)^(1/2))/C3 - (C4*Da*b*exp((2*(C3*w(m)*1i)^(1/2))/C3)*((2*(Da*w(m)*1i)^(1/2))/(Da*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))*(1 + (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*2i)/(w(m)*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))) - (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*1i)/w(m))) - (4*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))*(1 + (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*2i)/(w(m)*(exp((C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)) + exp(-(C1*Dox*w(m)*1i)^(1/2)/(C1*Dox)))) - (C2*Da*((Da*w(m)*1i)^(1/2)/Da - (2*(Da*w(m)*1i)^(1/2))/(Da*(exp((2*(Da*w(m)*1i)^(1/2))/Da) + 1)))*1i)/w(m))))*(C3*w(m)*1i)^(1/2))/C3)*1i)/(b*w(m)*(exp((2*(C3*w(m)*1i)^(1/2))/C3) + 1)) - ((eta0*w(m)*1i + eta0*w(m)*exp((2*(C3*w(m)*1i)^(1/2))/C3)*1i)*(C3*w(m)*1i)^(1/2)*1i)/(C3*b*w(m)*(exp((2*(C3*w(m)*1i)^(1/2))/C3) + 1)))*(exp((2*(C3*w(m)*1i)^(1/2))/C3) + 1)); 
end
for m=1:n
    Zreal(m)=(real(z(m)));
    Zimag(m)=(imag(z(m)));
end

plot(Zreal,-Zimag)
hold on
xlabel('Z_{real}')
ylabel('-Z_{imaginary}')
title('Nyquist plot for R_a = 10^{-5} cm and H_a \in [50,100] nm')
legend('H_a = 50nm','H_a = 60nm','H_a = 70nm','H_a = 80nm','H_a = 90nm','H_a = 100nm','location','northwest');
end