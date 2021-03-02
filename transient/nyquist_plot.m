figure
hold on
load Zt_100.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zt_200.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zt_500.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zt_1000.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zt_2000.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zt_5000.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zt_10000.mat
plot(-Zr(1:31),-Zj(1:31),'o')
load Zf
plot(-real(Z(1:31)),-imag(Z(1:31)),'ko')

