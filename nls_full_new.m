%memory cleanup:
clc;
clearvars;
close all;

%parameters:
L = 5; % optical fiber length
beta2 = -1; %dispersion coeffecient
gama=1; %nonlinear phase factor
NS = 1; % soliton order
nt = 1024; %number of FFT points
T = 32; %time axis window size
n = round(20*L*NS^2); % number of z steps
h = L/n; % step size in space (z)
dtau = (2*T)/nt; % step size in time axis 
tau = (-nt/2:nt/2-1)*dtau; % time array
omega = fftshift(-nt/2:nt/2-1)*(pi/T); % omega array
uin = sech(tau); %input pulse shape

%input pulse plot:
uin_f = fftshift(fft(uin)); %Fourier transform 
uin_pf = abs(uin_f).^2; % spectral power
uin_pf = uin_pf./max(uin_pf); %normalized spectral power
f = fftshift(omega)/(2*pi); %frequency array

%time domain plot:
subplot(2,1,1);
plot(tau, abs(uin).^2, '--ro');
hold on;
axis([-5 5 0 inf]);
xlabel('Normalized Time');
ylabel('Normalized Power');

%frequency domain plot:
subplot(2,1,2);
plot(f, uin_pf, '--ro');
hold on;
axis([-.5 .5 0 inf]);
xlabel('Normalized Frequency');
ylabel('Spectral Power');

%NLS operators:
D = exp(0.5i*beta2*omega.^2*h); %dispersion operator
N = gama*1i*(NS^2)*h; %nonlinear operator

ut = uin.*exp(abs(uin).^2.*N/2); %primary first half nonlinear step
for i=0:n
    uf = fft(ut).*D; %dispersion in Fourier domain
    ut = ifft(uf); %back to time domain
    ut = ut.*exp(abs(ut).^2.*N); %Nonlinear operator in time domain
end

%completing with last half step nonlinearity:
u_out = ut.*exp(-abs(ut).^2.*N/2); % Final field

%output pulse plot:
uf = fftshift(fft(u_out)); % Fourier transform
uout_pf = abs(uf).^2; % output spectrum power
uout_pf = uout_pf./max(uout_pf); % normalize output power

%time domain plot:
subplot(2,1,1)
plot(tau, abs(u_out).^2, '-b');
hold off;

%frequency domain plot:
subplot(2,1,2)
plot(f, uout_pf, '-b');
hold off;
