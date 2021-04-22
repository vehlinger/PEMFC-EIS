figure
load Zf_1V
plot(-real(Z),-imag(Z),'o')
hold on
load Zf_095V
plot(-real(Z),-imag(Z),'o')
load Zf_09V
plot(-real(Z),-imag(Z),'o')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')
legend('1.00 V', '0.95 V', '0.90 V','location','best')

figure
load Zf_075V
plot(-real(Z),-imag(Z),'o')
hold on
load Zf_072V
plot(-real(Z),-imag(Z),'o')
load Zf_07V
plot(-real(Z),-imag(Z),'o')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')
legend('0.75 V', '0.72 V', '0.70 V','location','best')

figure
load Zf_kappa10.mat
plot(-real(Z),-imag(Z),'o')
hold on
load Zf_kappa05.mat
plot(-real(Z),-imag(Z),'o')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')
legend('\kappa = 0.1 S/cm','\kappa = 0.05 S/cm','location','best')

figure
load Zf_a1000.mat
plot(-real(Z),-imag(Z),'o')
hold on
load Zf_a500.mat
plot(-real(Z),-imag(Z),'o')
xlabel('Z_r (\Omega cm^2)')
ylabel('-Z_j (\Omega cm^2)')
legend('a = 1000 cm^{-1}','a = 500 cm^{-1}','location','best')







