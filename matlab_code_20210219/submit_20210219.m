clc
clear all
close all

%% input field parameters
s = 2.5;                                             %aperture size unit:mm
pitch = 0.005;                                       % sampling interval
n0 = round(s/pitch);
r = s/pitch/2;
t = zeros(n0,1);
t (round(n0/2)-r+1:round(n0/2)+r,1) = 1;            %aperture
lam = 500e-6;                                       %illumination wavelength
k = 2*pi/lam;
n1 = size(t,1);
l = n1*pitch;
x1 = linspace(-l/2,l/2-pitch,n1)';                  %input plane coordinate
fx1 = linspace(-1/2/pitch,1/2/pitch-1/l,n1)';       %spatial frequency
figure,plot(x1,t);title('object')
t_FT = fftshift(fft(fftshift(t)));                  %Fourier spectrum of the aperture
figure,plot(fx1,abs(t_FT));title('FT of object')

%% output field parameters
    theta = 30*pi/180;                             %illumination angle
    z = 200;                                       %propagation distance
    
    S = s+0.02*z;                                  %size of the diffraction field
    L = (s+S)/2;
    x0 = z*tan(theta);                             %off-axis offset
    
    theta_min = asin(sin(theta)-lam/2/pitch);
    theta_max = asin(sin(theta)+lam/2/pitch);

    d1 = S/2-n0*pitch/2+x0-z*tan(theta_min);
    d2 = S/2-n0*pitch/2-x0+z*tan(theta_max);
    pitch_f = 1/(n0*pitch+max(d1,d2));            %sampling interval of the transfer function
    
    f_max = (x0+L)/lam/sqrt((x0+L)^2+z^2)-sin(theta)/lam;
    f_min = (x0-L)/lam/sqrt((x0-L)^2+z^2)-sin(theta)/lam;
    
    f_ul = min(f_max,1/2/pitch);
    f_ll = max(f_min,-1/2/pitch);
    M = ceil((f_ul-f_ll)/pitch_f);
    
    fx = f_ll+(1:M)*pitch_f-pitch_f;
    fx = fx';
    H = exp(1i*k*(lam*x0*fx+z*sqrt(1-(fx*lam+sin(theta)).^2))); % transfer function
    
    
    pitch_2 = round(1000/2/max(abs(f_ul),abs(f_ll)))/1000;      % sampling interval of the diffraction field
    n2 = ceil(S/pitch_2);
    S = n2*pitch_2;
    X = linspace(-S/2+x0,S/2+x0-pitch_2,n2)';
    K1 = n1/2/max(abs(fx1));
    rr = S/s;
    
%% NUFFT calculation
iflag = -1;
eps = 10^(-12);

tic
t_asmNUFT = nufft1d3(n1,x1/max(abs(x1))*pi,t,iflag,eps,M,(fx)*K1);
t_pro_asmNUFT = nufft1d3(M,(fx)*K1*rr,(H.*t_asmNUFT),-iflag,eps,n2,(X-x0)/(max(abs(X-x0)))*pi);
toc
t_pro_asmNUFT = t_pro_asmNUFT/max(abs(t_pro_asmNUFT))/sqrt(1i);

figure,plot(X-x0,abs(t_pro_asmNUFT));set(gca,'looseInset',[0 0 0 0]);title('proposed amplitude')
figure,plot(X-x0,angle(t_pro_asmNUFT));set(gca,'looseInset',[0 0 0 0]);title('proposed phase')

%% analytical integral 
uu = zeros(n2,1);
tic
for j = 1:n2
      fun = @(xn) 1/2/pi*z./sqrt((X(j)-xn).^2+z^2).*(1./sqrt((X(j)-xn).^2+z^2)...
               -1i*k/pi).*exp(1i*k*sqrt((X(j)-xn).^2+z^2))./sqrt((X(j)-xn).^2+z^2).*exp(1i*k*xn*sin(theta));
    uu(j,1) = integral(fun,-(r-0)*pitch,(r-1)*pitch);
end
toc
uu = uu/max(abs(uu))./exp(1i*k*X*sin(theta));

phase_rsi = (angle(uu));
amplitude_rsi = abs(uu);
intensity_rsi = amplitude_rsi.^2;

figure,plot(X-x0,(phase_rsi));title('Analytical integral phase ')
set(gca,'looseInset',[0 0 0 0])
figure,plot(X-x0,amplitude_rsi);title('Analytical integral amplitude')
set(gca,'looseInset',[0 0 0 0])

%% SNR
aa = (t_pro_asmNUFT);
bb = (uu);
alpha_be = (sum(aa.*conj(bb)))/(sum((abs(bb)).^2));
snr_be = 10*log10((sum((abs(aa)).^2))/(sum((abs(aa-alpha_be*bb)).^2)))



