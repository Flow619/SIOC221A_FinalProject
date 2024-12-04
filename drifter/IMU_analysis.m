% for tracking IMU drift 

%load in data to table; if multiple files concatenate vertically

%data1 = readtable('15X09X48.CSV');
% tank1 = readtable('data_dep1/13X11X12');
% tank2 = readtable('07X15X49');
% tank3 = readtable('07X15X52');
clear

data2 = readtable('14X11X21.CSV');
data3 = readtable('14X11X24.CSV');
data4 = readtable('14X11X27.CSV');
data5 = readtable('14X11X30.CSV');
data6 = readtable('14X11X33.CSV');
data7 = readtable('14X11X34.CSV');
data8 = readtable('14X11X36.CSV');
data9 = readtable('14X11X39.CSV');
data10 = readtable('14X11X42.CSV');
data11 = readtable('14X11X45.CSV');
data12 = readtable('14X11X48.CSV');
data13 = readtable('15X11X06.CSV');
data14 = readtable('15X11X09.CSV');
data15 = readtable('15X11X12.CSV');
data16 = readtable('15X11X15.CSV');
data17 = readtable('15X11X18.CSV');
data18 = readtable('15X11X21.CSV');
data19 = readtable('15X11X24.CSV');
data20 = readtable('15X11X27.CSV');
data21 = readtable('15X11X30.CSV');
data22 = readtable('15X11X34.CSV');
data23 = readtable('15X11X37.CSV');
data24 = readtable('15X11X40.CSV');
datatable = [ data2; data3; data4; data5; data6; data7; data8; data9; data10; data11; data12; 
    data13; data14; data15; data16; data17; data18; data19; data20; data21; data22; data23; data24];
save('Deploy2.mat', 'datatable');
%%
% figure(1)
% clf
% tank3time = datetime([2024, 03, 07, 15, 52, 00]);
% tank2time = datetime([2024, 03, 07, 15, 49, 00]);
tank1time = datetime([2024, 12, 04, 10, 37, 00]);
acc1 = table2array(datatable(:,4:6)).*9.81./1000;
% acc2 = table2array(tank2(100:900,4:6)).*9.81./1000;
% acc3 = table2array(tank3(1:250,4:6)).*9.81./1000;
% acc1_1 = acc1(200:750,:);
% acc1_2 = acc(800:end-150,:);
time1 = table2array(datatable(:,1))./1000;
% time2 = table2array(tank2(100:900,1))./1000;
% time3 = table2array(tank3(1:250,1))./1000;
% tank2timeaxis = datetime((time2/60/60/24 + datenum(tank2time)), 'convertfrom','datenum');
% tank3timeaxis = datetime((time3/60/60/24 + datenum(tank3time)), 'convertfrom','datenum');
tank1timeaxis = datetime((time1/60/60/24 + datenum(tank1time)), 'convertfrom','datenum');
tanktime1 = tank1timeaxis(200:750);
% tanktime2 = tank1timeaxis(800:end-150);

%%
figure(1);clf
subplot(211)
plot(tank1timeaxis,detrend(acc1))
title('Test 1 IMU Acceleration Timeseries')
ylabel('Acceleration [m/s^2]')
%%

%
%sampling_rate = nanmean(diff(tanktime1));
dt = 1/8;%sampling_rate;
Fs = 1/dt;
nbands = 6;
alpha = 0.1;
%
subplot(212)
%take the fourier transform of the data
[Sx,Su,Sl,f,fsel] = autospectrum(acc1(:,1), dt, nbands, alpha);
loglog(f,Sx)
grid on; hold on
[Sy,Su,Sl,f,fsel] = autospectrum(acc1(:,2), dt, nbands, alpha);
loglog(f,Sy)
grid on; hold on
[Sz,Su,Sl,f,fsel] = autospectrum(acc1(:,3), dt, nbands, alpha);
figure(2);
loglog(f,Sx);
hold on
loglog(f,Sz);
grid on; hold on
loglog([1 1], ylim ,'-b');
title('Wave Tank Test 1 Acceleration Spectra')
xlabel('Frequency [Hz]')
ylabel('\Phi_a [(m/s^2)^2/Hz]')
%%

data2 = readtable('13X11X06.CSV');
data3 = readtable('13X11X09.CSV');
data4 = readtable('13X11X12.CSV');
data5 = readtable('13X11X15.CSV');
data6 = readtable('13X11X18.CSV');
data7 = readtable('13X11X21.CSV');
data8 = readtable('13X11X24.CSV');
data9 = readtable('13X11X27.CSV');
data10 = readtable('13X11X30.CSV');
data11 = readtable('13X11X33.CSV');
data12 = readtable('13X11X36.CSV');
data13 = readtable('13X11X39.CSV');
data14 = readtable('13X11X42.CSV');
% data15 = readtable('13X11X06.CSV');
% data16 = readtable('13X11X06.CSV');

datatable = [ data2; data3; data4; data5; data6; data7; data8; data9; data10; data11; data12; data13; data14;];
% clear data1 data2 data3 data4 data5 data6 data7 data8 data9 data10 data11 data12 data13 data14
%% Look at the data
%start = 700;
%ennd = 1050;
start = 2000;
ennd = size(datatable,1);

%data(6400:6700,:) = 0;

%accelerometer in g
%gyro in degrees / s
%
starttime = datetime([2024, 03, 15, 10, 51, 00]);
acceleration_data = detrend(table2array(datatable(start:ennd, 4:6)))*9.81./1000;  % Assuming columns 1, 2, and 3 correspond to X, Y, and Z axes
time = table2array( datatable(start:ennd,1))./1000;  %-  table2array(data(start, 1))/1000;
timeaxis = datetime((time/60/60/24 + datenum(starttime)), 'convertfrom','datenum');
timeaxis1 = timeaxis(1:4350); timeaxis2 = timeaxis(4750:end);
%acceleration_data(:,3) = acceleration_data(:,3) %- mean(acceleration_data(:,3));
gyro_data = table2array(datatable(start:ennd, 7:9));
gyro1 = gyro_data(1:4350,:); 
gyro2 = gyro_data(4750:end,:);
compass_data = detrend(table2array(datatable(start:ennd, 10:12)));

acc_data1 = acceleration_data(1:4350,:); time1 = time(1:4350);
acc_data2 = acceleration_data(4750:end,:); time2 = time(4750:end);

figure(1); clf
subplot(3,1,1)
plot(timeaxis, acceleration_data);
title('IMU Acceleration Timeseries');
ylabel('Acceleration (m/s^2)');
legend('X-axis', 'Y-axis', 'Z-axis');
grid on;

subplot(3,1,2)

plot(timeaxis, gyro_data);
title('IMU Gyroscope Timeseries')
ylabel('Gryo Degrees/s');
legend('x', 'y', 'z')
ylim([-100, 100])
subplot(3,1,3)

plot(timeaxis, detrend(compass_data));
title('IMU Compass Timeseries')
ylabel('Compass Degrees')
legend('x' , 'y', 'z')
xlabel('Time')


%% Highpass filter
dt = nanmean(diff(time1));
%Fc = 0.001;%frequency cutoff between 0 and 0.5 (normalized freq with fs = 1)
N = 3;% number of coeffs
n = -(N-1)/2:(N-1)/2 + 1;
window  = hamming(length(acc_data1));
fc = 1/300;
accf = highpass(acc_data1, fc, (1/dt), 'steepness', 0.9);
%hp = lowpass(hp, 0.05, 1/dt);

%impulse response of the filter
figure(2)
clf
plot(time1, detrend(accf))
grid on
title('highpass filter, fc = 1/300')
% highpass filter using butterworth
%
fs = 5;
order = 6;
[b, a] = butter(order, 0.2, 'low');
accf = filtfilt(b, a, acc_data1);

figure(3)
clf
plot(time1, accf);
title('butter filter, fc = 1/300')
% 
% 
% highpass filter plotted
[h1,h1_windowed,n] = idl_highPass(pi*0.4,20);

w = (0:0.01:pi);
H1 = freqz(h1,1,w);
H1_windowed = freqz(h1_windowed,1,w);

Hh = hamming(length(acc_data1));

figure(4)
clf
plot(w/pi,abs(H1));
hold on
grid on
plot(w/pi,abs(H1_windowed));
xlabel('Normalized frequency ( 1 = Nyquist (fs/2) )')
legend('rectangular windowed', 'hamming window')
title('high pass filter')

%%
accf = acc_data1;
t_mod = time1;
wn = 12;
N = length(accf);
while(mod(N,wn)~=0)
    N = N-1;
    t_mod = t_mod(1:N);
    accf = accf(1:N,:);
end
%

data = accf(:,1);
chunk = length(accf)/wn;
[fx,ax,Parsevalx] = spectrumCB(t_mod, data, chunk);
data = accf(:,2);
chunk = length(accf)/wn;
[fy,ay,Parsevaly] = spectrumCB(t_mod, data, chunk);
data = accf(:,3);
chunk = length(accf)/wn;
[fz,az,Parsevalz] = spectrumCB(t_mod, data, chunk);



figure(6)
clf
loglog(fx, ax);
hold on
loglog(fy, ay);
hold on
loglog(fz,az);
grid on
hold on
f_5s = 1/5.88;
f_11s = 1/11.5;
loglog([f_5s f_5s], ylim ,'-b');
loglog([f_11s f_11s], ylim ,'-b');
text(f_5s, 30000, '5 second wave')
text(f_11s, 50000, '11 second wave')
title('Scripps Pier Drifter Acceleration Spectra')
xlabel('Frequency [Hz]')
ylabel('\Phi_a [(m/s^2)^2/Hz]')
%% get velocity

decim = 1;
Fs = 1/5;
fuse = imufilter('SampleRate',Fs*4,'DecimationFactor',decim);
[Ornt, Av] = step(fuse, acc_data1, gyro1);
quat = fuse(acc_data1,gyro1);
timed = (0:decim:size(acc_data1,1)-1)/Fs;

Oq = eulerd(Ornt,'ZYX','frame');

Corrected_acc2 = Oq.*acc_data1*Oq';
%%
figure(4);clf
plot(timeaxis1,Corrected_acc2)
title('Orientation Estimate')
legend('Z-axis', 'Y-axis', 'X-axis')
xlabel('Time (s)')
ylabel('Rotation (degrees)')
%%
fc = 0.005;
dt = 0.2;
fs = 1/dt;
[B,A] = butter(1,2*fc/fs, 'high');

figure(4);clf
disp = cumtrapz(detrend(Av));
plot(timeaxis1, disp);
title('Position Estimate')
legend('Z-axis', 'Y-axis', 'X-axis')
xlabel('Time (s)')
ylabel('Position (m)')


%%

%%

velocity = detrend(cumtrapz(t_mod, accf));

disp = detrend(cumtrapz(t_mod, velocity));

figure(4)
clf
plot(t_mod, detrend(velocity(:,1)))

%%
find_nearest = @(array, value) array(abs(array - value) == min(abs(array - value)));
%make sure accelerometer is taking burst samples 0.25 second sampling rate over a 60 second burst
% make sure sampling is evenly spaced -- good clock
%make sure data is telemetered
%fourier analysis:
%   x-axis: frequency
%   y-axis: velocity squared / frequency

% use fast fourier transform (fft) 
% need timeseries data, bins of frequency bands (1hz - 0.05hz for ocean waves)
% nyquist frequency of 1hz,
% minimum frequency of 2hz or more at inner bound
% minimum frequency of 0.05 hz at lower outer band

sampling_rate = nanmean(diff(time1));
dt = 0.2;%sampling_rate;
Fs = 1/dt;
nbands = 20;
alpha = 0.1;
%
figure(3); clf
%take the fourier transform of the data
[Sx,Su,Sl,f,fsel] = autospectrum(acc_data2(:,1), dt, nbands, alpha);
loglog(f,Sx)
%hold on;
%loglog(Df,DSx)
xlabel('Frequency [Hz]','FontSize',14)
ylabel('Acceleration Spectrum [g^2/Hz]','Fontsize',14)
grid on;
hold on
%[DSy,DSu,DSl,Df,Dfsel] = autospectrum(dataTable.RawY, Dt, nbands, alpha);   %
[Sy,Su,Sl,f,fsel] = autospectrum(acc_data2(:,2), dt, nbands, alpha);
loglog(f,Sy)
hold on;
%loglog(Df,DSy)
%hold on

%[DSz,DSu,DSl,Df,Dfsel] = autospectrum(dataTable.RawZ, Dt, nbands, alpha);   %
[Sz,Su,Sl,f,fsel] = autospectrum(acc_data2(:,3), dt, nbands, alpha);
loglog(f,Sz)
%hold on;
%loglog(Df,DSz)

legend('Acc_X', 'Acc_Y', 'Acc_Z', 'location', 'best');

%%
Ax = acc_data1(:,1);
Ay = acc_data1(:,2);
Az = acc_data1(:,3);
fAx = fftshift(fft(Ax));
fAy = fftshift(fft(Ay));
fAz = fftshift(fft(Az));

N = size(acc_data1,1);
dt = 1/5;
fN = dt/2;
fn = 1/(2*dt);
df = 1/(N*dt);
f = (-fn:df:fn-df)';
f = f(1:N);

kig = find(abs(f) < 1/25 & abs(f) > 1/900);
kss = find(abs(f) < 1/20);


Axss = fAz;
Axss(kss) = 0;
Axssreconstruct = real(ifft(ifftshift(Axss)));
Hxss = 4*std(Axssreconstruct) / 9.81;


figure(4)
clf
loglog(f, X)


%%

Axig = fAx;
Axig(kig) = 0;
Aigreconstruct = real(ifft(ifftshift(Axig)));
Higx = 4*std(Aigreconstruct) / 9.81;
fAy(kig) = 0;
fAz(kig) = 0;

fAx2 = ifftshift(fAx);
fAy2 = ifftshift(fAy);
fAz2 = ifftshift(fAz);

Axs = real(ifft(fAx2));
Ays = real(ifft(fAy2));
Azs = real(ifft(fAz2));
figure(4);
clf
subplot(311)
plot(time, Ax, 'b', time, Ax-Axs);
subplot(312)
plot(time, Ay, 'b', time, Ay-Ays);
subplot(313)
plot(time, Az, 'b', time, Az-Azs);


figure(5);
clf
[Sx,Su,Sl,f,fsel] = autospectrum(Ax-Axs, dt, nbands, alpha);
loglog(f,Sx)
hold on
[Sy,Su,Sl,f,fsel] = autospectrum(Ay-Ays, dt, nbands, alpha);
loglog(f,Sy)
hold on
[Sz,Su,Sl,f,fsel] = autospectrum(Az-Azs, dt, nbands, alpha);
loglog(f,Sz)
ylim([1 100000])

%%

% integrate the displacement over the frequency f = 11s
N = size(acceleration_data,1);
Accz = acceleration_data(:,3);

%filter data
fN = dt/2;
fn = 1/(2*dt);
df = 1/(N*dt);
f = (-fn:df:fn-df)';
f = f(1:N);
kwave = find(abs(f) > 1/12 & abs(f) < 1/9);
Acczf = fftshift(fft(Accz));

velocity = cumtrapz(Acczf(kwave));
disp = cumtrapz(velocity);

figure(4)
clf
plot( disp);


%%
function [hd,hd_windowed,n] = idl_highPass(wc,N)

if mod(N,2)==0
  a = N/2; 
  delta_n = [zeros(1,a) 1 zeros(1,a-1)];
else
  a = (N-1)/2;
  delta_n = [zeros(1,a) 1 zeros(1,a)];
end

m = 0:1:N-1;

n = m-a+eps;

hd = delta_n - sin(wc.*n)./(pi.*n);

% you should window this:
hd_windowed = hd .* hamming(N)';

end
%%

% needs to aqcuire data in bursts over a certain frequency 
function [S,Su,Sl,f,fsel] = autospectrum(x,dt,nbands,alpha)
%Computes the autospectrum of input time series x
%Input:
%x: time series
%dt: sample period
%nbands: number of frequency bands to smooth
%alpha: (1-alpha)*100 confidence intervals
%Output:
%S: autospectrum
%Su: 95% confidence interval upper bound
%Sl: 95% confidence interval lower bound
%f: frequency, units = cycles/dt
%ensure that length of x is an even number
N = length(x);
if(mod(N,2)~=0)
N = N-1;
x = x(1:N);
end
%Fourier transform
X = fft(x);
%one-sided frequencies, ignore the mean
f = (1:N/2)'/(N*dt);
fsel = find(abs(f)>1/000);
%one-sided spectral density
S = (abs(X(2:N/2+1)).^2)*dt/N;
%double amplitudes except Nyquist
S(1:end-1) = 2*S(1:end-1);
%smooth with a running mean
w = ones(nbands,1)/nbands;
%S = conv(S,w,'valid');
%f = conv(f,w,'valid');
%confident intervals
nu = 2*nbands;
Su = S*nu/chi2inv(alpha/2,nu);
Sl = S*nu/chi2inv(1-alpha/2,nu);

%inverse fft
%I = ifft(x)


end

%%
function [f,X,Parseval] = spectrumCB(t_mod, data, chunk)

ind1 = 1;

for i=1:floor(chunk\length(data))*2
    
    if ind1+chunk<length(data)
        data1(:,i) = data(ind1:ind1+chunk).*hamming(chunk+1);
    end
    ind1 = ind1+chunk/2;
end
%frequencies
dt = nanmean(diff(t_mod)); %time between samples
fn = 1/2/dt; %Nyquist Frequency
N = length(data1); 
T = dt*N; %time length of record
df = 1/T; %fundamental frequency
f = 0:df:fn; %frequency vector


for i=1:length(data1(1,:))
    data2(:,i) = detrend(data1(:,i));
    X = fft(data2(:,i));
    amp = abs(X(1:(N+1)/2)).^2;
    amp = amp / N.^2;
    amp = amp .*2;
    amp = amp / df;
    A(:,i) = amp;
end
%average fourier transform chunks
X = mean(A,2);
variance = var(data);
int_spec = trapz(f, X);
Parseval = int_spec/variance;

end