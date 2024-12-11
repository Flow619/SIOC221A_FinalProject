close all;
clear;
clc;

% Author: Gabriel Gekas

deployment = 1; % choose deployment

dirpath = strrep(pwd, 'Gabe_Data_Analysis', 'drifter'); % get data directory
dirs = {'/data_dep1', '/data_dep2'}; % define deployments
DIR = [dirpath, dirs{deployment}];

csvFiles = dir(fullfile(DIR, '*.CSV')); % gets all .csv files
for i = 1:length(csvFiles)
    % get input file :
    filepath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    data{i} = readtable(filepath);
end

start_index = [3, 13]; % start of data file deployment
end_index = [11, length(csvFiles)-1]; % end of data file deployment

fe = 5; % sampling rate in Hz
fcut = 1/40; % cutoff frequencies for highpass filter
order = 4; % order of butterworth highpass filter
count = 1;


figure; hold on;
for i = start_index(deployment):end_index(deployment)
    plot((i-1)*length(data{i}.accZ) + [1:length(data{i}.accZ)], data{i}.accZ)

    ACC = [data{i}.accX, data{i}.accY, data{i}.accZ]';
    [acc{count}, vel{count}, disp{count}] = compute_DISP_high_pass(ACC, fe, fcut, order);

    count = count + 1;
end



% %% plot the filter :
% Ordercut=order;
% Wcut=2*fcut/fe;
% [b,a] = butter(Ordercut,Wcut,'high');
% freqz(b,a,[],fe)


dt = 1/fe;
nseg = 5;

%% acceleration :
for i = 1:count-1
    vec = acc{i}';
    x = sqrt(sum(vec.^2, 2));
    [Pout{i}, fout{i}, errbars{i}] = gabespectra(x, dt, nseg);
end

figure; a = tiledlayout(2, 1);
nexttile(1); hold on;
nexttile(2); hold on;
for i = 1:count-1
    nexttile(1);
    plot(fout{i}, Pout{i})
    nexttile(2);
    plot(fout{i}, Pout{i})

    if i == 1
        Pmat = Pout{i};
    else
        Pmat = [Pmat, Pout{i}];
    end
end
nexttile(1);
mean_line = plot(fout{1}, mean(Pmat, 2), 'k-', 'LineWidth', 1.4);
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
legend([mean_line], {'Mean FFT'})

nexttile(2);
mean_line = plot(fout{1}, mean(Pmat, 2), 'k-', 'LineWidth', 1.4);
set(gca, 'XScale','log', 'YScale','log')
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
legend([mean_line], {'Mean FFT'})
title(a, ['$\sqrt{accX^2 + accY^2 + accZ^2}, \:', num2str(nseg), '\: segments', ', \: Deployment: \:', num2str(deployment), '$'], 'Interpreter', 'latex')


%% velocity :
for i = 1:count-1
    vec = vel{i}';
    x = sqrt(sum(vec.^2, 2));
    [Pout{i}, fout{i}, errbars{i}] = gabespectra(x, dt, nseg);
end

figure; a = tiledlayout(2, 1);
nexttile(1); hold on;
nexttile(2); hold on;
for i = 1:count-1
    nexttile(1);
    plot(fout{i}, Pout{i})
    nexttile(2);
    plot(fout{i}, Pout{i})

    if i == 1
        Pmat = Pout{i};
    else
        Pmat = [Pmat, Pout{i}];
    end
end
nexttile(1);
mean_line = plot(fout{1}, mean(Pmat, 2), 'k-', 'LineWidth', 1.4);
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
legend([mean_line], {'Mean FFT'})

nexttile(2);
mean_line = plot(fout{1}, mean(Pmat, 2), 'k-', 'LineWidth', 1.4);
set(gca, 'XScale','log', 'YScale','log')
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
legend([mean_line], {'Mean FFT'})
title(a, ['$\sqrt{velX^2 + velY^2 + velZ^2}, \:', num2str(nseg), '\: segments', ', \: Deployment: \:', num2str(deployment), '$'], 'Interpreter', 'latex')


%% displacement :

for i = 1:count-1
    vec = disp{i}';
    x = sqrt(sum(vec.^2, 2));
    [Pout{i}, fout{i}, errbars{i}] = gabespectra(x, dt, nseg);
end

figure; a = tiledlayout(2, 1);
nexttile(1); hold on;
nexttile(2); hold on;
for i = 1:count-1
    nexttile(1);
    plot(fout{i}, Pout{i})
    nexttile(2);
    plot(fout{i}, Pout{i})

    if i == 1
        Pmat = Pout{i};
    else
        Pmat = [Pmat, Pout{i}];
    end
end
nexttile(1);
mean_line = plot(fout{1}, mean(Pmat, 2), 'k-', 'LineWidth', 1.4);
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
legend([mean_line], {'Mean FFT'})

nexttile(2);
mean_line = plot(fout{1}, mean(Pmat, 2), 'k-', 'LineWidth', 1.4);
line1 = plot(fout{1}, fout{1}.^(-4), 'r-', 'LineWidth', 1.4);
% line2 = plot(fout{1}, fout{1}.^(-3), 'b-', 'LineWidth', 1.4);
set(gca, 'XScale','log', 'YScale','log')
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
legend([mean_line, line1], {'Mean FFT', 'f^{-4}'})
title(a, ['$\sqrt{X^2 + Y^2 + Z^2}, \:', num2str(nseg), '\: segments', ', \: Deployment: \:', num2str(deployment), '$'], 'Interpreter', 'latex')


function [Pout, fout, errbars] = gabespectra(x, dt, nseg)
% computes filtered power spectrum of signal when given the signal,
% timestep, and number of windows to filter data through.

    % taken from Jennifer Mackinnon :
    if nargin < 3; nseg = 10; end % if you choose not to input a number of segments
    N = length(x); % length of original time series
    Nseg = floor(N/nseg); % how many segments will fit
    xshort = x(1:nseg*Nseg); % slightly shorter time series of a size that is an integer
    xreshape1 = reshape(xshort, Nseg, nseg); % sequential segments
    xreshape2 = reshape(xshort((floor(Nseg/2)+1):(end-ceil(Nseg/2))), Nseg, nseg-1); % overlapping segments between sequential segments
    xreshape = [xreshape1, xreshape2]; % combine all the segments. (Nseg x (2*nseg-1))
    xdet = detrend(xreshape); % detrends in one go
    wind = hann(Nseg)*ones(1, 2*nseg-1); % hanning windows
    xwindowed = xdet.*wind; % multiplies by hanning windows
    xf = fftshift(fft(xwindowed)); % ffts each column


    % written for the assignment :
    df = 1 / (Nseg*dt); % fundamental frequency for each segment
    f = ceil(-Nseg/2:Nseg/2-1)'*df; % get frequency vector for each segment
    P = (dt/Nseg) * abs(xf).^2; % square the amplitudes and normalize
    % flip and add coeffs from negative frequencies to positive frequencies
    if mod(Nseg, 2) == 0 % even
        P2 = flipud(P(1:ceil(Nseg/2)+1, :));
        P2(2:(end-1), :) = P2(2:(end-1), :) + P((ceil(Nseg/2)+2):end, :);
    else % odd
        P2 = flipud(P(1:ceil(Nseg/2), :));
        P2(2:(end-1), :) = P2(2:(end-1), :) + P((ceil(Nseg/2)+1):end-1, :);
    end
    fout = flipud(abs(f(f <= 0))); % get frequencies 0:fn

    Pout = mean(P2, 2); % get power spectra average across columns
    Pout = Pout / mean(hann(Nseg).^2); % adjust for hanning window power reduction

    % calculate error bars
    nu = nseg;
    errhigh = nu / (chi2inv(0.025, nu));
    errlow = nu / (chi2inv(0.975, nu));
    errbars = [errlow, errhigh];
    
end