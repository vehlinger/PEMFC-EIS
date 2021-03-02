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




