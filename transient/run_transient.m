clear all; close all; clc

n = 6;               % number of unknowns at each mesh point
nj = 21;             % number of mesh points
C = zeros(nj,n);     % change variable 

% parameters
alpha = 1;  
sigma = 7;   % S/cm  
kappa = 0.1; % S/cm
a = 1e3;     % 1/cm  
Cdl = 0.2;   % F/cm2 
D = 0.3;     % cm2/s (O2 in water)
params = [n nj alpha sigma kappa a Cdl D];

% operating conditions
L = 0.001;   % cm
T0 = 353.15; % K
Vcell = 1; % V
RH = 0.5;
p = 1;
Pwsat = exp(11.6832-3816.44/(T0-46.13)); 
C0 = 0.21*(1-RH*(Pwsat/p));   
op_cond = [L T0 Vcell C0 p];

load C_ss.mat C_ss
C_ss = steady_state(C_ss,n,nj,params,op_cond);

%% transient 
frange = logspace(-3,-2,11); % method goes unstable around 60 mHz
for ii = 1:length(frange)
    f = frange(ii); % Hz
    disp(f)
    T = 1/f;        % period (s)
    omega = 2*pi*f; % angular frequency
    deltaV = 1e-5;  % V
    n_cycle = 5;

    tfinal = T*n_cycle;
    dt = 0.0001*T;

    time(1) = 0;
    k = 1;

    C = C_ss;
    Cp = C_ss;
    Ct(:,:,1) = reshape(C_ss,[nj,n,1]);

    while time < tfinal
        k = k+1;
        time(k) = time(k-1)+dt;
        op_cond(3) = Vcell+deltaV*cos(omega*time(k));
        C = transient(C,n,nj,params,op_cond,Cp,dt,time(k));
        Ct(:,:,k) = C;
        Cp = C;
    end

    V = reshape(Ct(1,2,:),1,size(Ct,3));
    i = reshape(Ct(end,1,:),1,size(Ct,3));
    
    figure(1)
    yyaxis left
    plot(time,V)
    yyaxis right
    plot(time,i)
    
    Ir(ii) = trapz(i.*cos(omega.*time))/tfinal;
    Ij(ii) = -trapz(i.*sin(omega.*time))/tfinal;
    Vr(ii) = trapz(V.*cos(omega.*time))/tfinal;
    Vj(ii) = -trapz(V.*sin(omega.*time))/tfinal;
    
    Zr(ii) = real((Vr(ii)+1j*Vj(ii))/(Ir(ii)+1j*Ij(ii)));
    Zj(ii) = imag((Vr(ii)+1j*Vj(ii))/(Ir(ii)+1j*Ij(ii)));
    
    clear time V i C Cp Ct
end 

%%
% note Zr must be mulitplied by -1 to get positive real impedance
figure
subplot(1,2,1)
plot(-Zr,-Zj,'o-')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')

subplot(1,2,2)
yyaxis left
loglog(frange,-Zr)
ylabel('Z_r (\Omega cm^2)')
yyaxis right
loglog(frange,-Zj)
xlabel('frequency (Hz)')
ylabel('-Z_j (\Omega cm^2)')
