clear all; close all; clc

timerVal = tic;
n = 6;               % number of unknowns at each mesh point 
nj = 21;             % number of mesh points
C = zeros(nj,n);     % change variable 

% parameters
alpha = 1;
sigma = 7;   % S/cm
kappa = 0.05; % S/cm
a = 1000;     % 1/cm 
Cdl = 0.2;   % F/cm2 
D = 0.3;     % cm2/s (O2 in water)
params = [n nj alpha sigma kappa a Cdl D];

% operating conditions
L = 0.001;   % cm
T0 = 353.15; % K
Vcell = 1; % V
deltaV = 1e-5;
RH = 0.5;
p = 1; % bar
Pwsat = exp(11.6832-3816.44/(T0-46.13)); 
C0 = 0.21*(1-RH*(Pwsat/p)); % mol/cm3  

for j = 1:length(Vcell)
    op_cond = [L T0 Vcell(j) C0 p];

    load C_ss.mat C_ss
    C_ss = steady_state(C_ss,n,nj,params,op_cond);

    frange = logspace(-3,6,91);

    for ii = 1:length(frange)
        f = frange(ii);       % frequency (Hz)
        omega = 2*pi*f; % angular frequency
        op_cond = [L T0 deltaV C0 p omega];
        Ctilde = complex(C_ss);
        Ctilde = freq_response(Ctilde,n,nj,params,op_cond,C_ss);
        Z(ii) = Ctilde(end,2)/(Ctilde(end,1));
    end
    filename = sprintf('Z%umV.mat',Vcell(j)*1000);
    save(filename, 'frange','Z')
end

tend = toc(timerVal);

%%
% note Zr must be mulitplied by -1 to get positive real impedance
figure
subplot(1,2,1)
plot(-real(Z),-imag(Z),'o-')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')

subplot(1,2,2)
yyaxis left
loglog(frange,-real(Z))
ylabel('Z_r (\Omega cm^2)')
yyaxis right
loglog(frange,-imag(Z))
xlabel('frequency (Hz)')
ylabel('-Z_j (\Omega cm^2)')
