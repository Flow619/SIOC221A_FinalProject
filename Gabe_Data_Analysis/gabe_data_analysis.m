close all;
clear;
clc;

% Author: Gabriel Gekas

deployment = 1; % choose deployment

dirpath = strrep(pwd, 'Gabe_Data_Analysis', 'drifter'); % get data directory
dirs = {'/data_dep1', '/data_dep2'}; % define deployments
DIR = [dirpath, dirs{deployment}];
start_index = [3, 1]; % start of data file deployment
end_index = [11, 10]; % end of data file deployment

csvFiles = dir(fullfile(DIR, '*.CSV')); % gets all .csv files
for i = 1:length(csvFiles)
    % get input file :
    filepath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    data{i} = readtable(filepath);
end

count = 1;
for i = start_index(deployment):end_index(deployment)
    acc_x{count} = data{i}.accX;
    acc_y{count} = data{i}.accY;
    acc_z{count} = data{i}.accZ;
    count = count + 1;
end

dt = 0.2;
nseg = 9;

for i = 1:count-1
    vec = [acc_x{i}, acc_y{i}, acc_z{i}];
    x = sqrt(sum(vec.^2, 2));
    [Pout{i}, fout{i}, errbars{i}] = gabespectra(x, dt, nseg);
end

figure; a = tiledlayout(2, 1);
nexttile; hold on;
for i = 1:count-1
    plot(fout{i}, Pout{i})
end
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');

nexttile; hold on;
for i = 1:count-1
    line_obj(i) = plot(fout{i}, Pout{i});
    % errneg = Pout{i}(end) - errbars{i}(1)*Pout{i}(end);
    % errpos = errbars{i}(2)*Pout{i}(end) - Pout{i}(end);
    % errorbar(fout{i}(end)+i, Pout{i}(end), errneg, errpos)
end
set(gca, 'XScale','log', 'YScale','log')
title('One-sided Power Spectrum');
ylabel('Power [$$ \frac{units^2}{Hz} $$]', 'Interpreter', 'latex');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');

title(a, '$\sqrt{accX^2 + accY^2 + accZ^2}, \: 9 \: segments$', 'Interpreter', 'latex')

% % satisfies Parseval's Theorem: these should be equal
% disp(['Energy in time domain : ', num2str(mean(x.^2))])
% for i = 1:length(nseg)
%     df = max(diff(fout{i}));
%     disp(['Windows : ', num2str(nseg(i))])
%     disp(['gabespectra(.) energy in frequency domain : ', num2str(sum(Pout{i}*df))])
% end



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