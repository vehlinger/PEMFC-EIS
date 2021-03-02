clear all; close all; clc

timerVal = tic;
n = 6;               % number of unknowns at each mesh point 
nj = 21;             % number of mesh points
C = zeros(nj,n);     % change variable 

% parameters
alpha = 1;
sigma = 7;   % S/cm
kappa = 0.1; % S/cm
a = 1e3;     % 1/cm 
Cdl = 2e-5;  % F/cm2 
D = 0.3;     % cm2/s (O2 in water)
params = [n nj alpha sigma kappa a Cdl D];

% operating conditions
L = 0.001;   % cm
T0 = 353.15; % K
Vcell = 0.72; % V
deltaV = 1e-5;
RH = 0.5;
p = 1;
Pwsat = exp(11.6832-3816.44/(T0-46.13)); 
C0 = 0.21*(1-RH*(Pwsat/p));   
op_cond = [L T0 Vcell C0 p];

load C_ss.mat C_ss
C_ss = steady_state(C_ss,n,nj,params,op_cond);

frange = logspace(-3,10,131);

for ii = 1:length(frange)
    f = frange(ii);       % frequency (Hz)
    omega = 2*pi*f; % angular frequency
    op_cond = [L T0 deltaV C0 p omega];
    Ctilde = complex(C_ss);
    Ctilde = freq_response(Ctilde,n,nj,params,op_cond,C_ss);
    Z(ii) = Ctilde(end,2)/(Ctilde(end,1));
end

tend = toc(timerVal);

figure
plot(-real(Z),-imag(Z),'o-')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')
