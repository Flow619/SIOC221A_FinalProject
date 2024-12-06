close all
clear all
clc

%% Load Data
load("FFT_Deploy1_Chunking.mat")
load("FFT_Deploy2_Chunking.mat")

%% Plot All FFTs
Colors = jet(57);

figure
hold on

% loop through and plot all FFTs
for i = 1:27
    plot(freq,Deploy1(:,i),"Color",[0,0,1,0.2],'HandleVisibility','off')
end

for i = 1:30
    plot(freq,Deploy2(:,i),"Color",[0,0,1,0.2],'HandleVisibility','off')
end

yscale('log')
xscale('log')

% Calculate Means
Mean_Deploy1 = mean(Deploy1,2);
Mean_Deploy2 = mean(Deploy2,2);
Mean_Deployments = mean([Deploy1,Deploy2],2);

% Plot Means
plot(freq,Mean_Deploy1,'r','LineWidth',2,'DisplayName','Mean Deploy #1')
plot(freq,Mean_Deploy2,'g','LineWidth',2,'DisplayName','Mean Deploy #2')
plot(freq,Mean_Deployments,'k-','LineWidth',2.5,'DisplayName','Mean Deploy #1 & #2')

set(gca,'FontSize',20)

l = legend
set(l,'fontsize',14)

title('Deployment #1 & #2 with Chunking')

%% Mean Deploy 1 & 2 FFT
figure
plot(freq,Mean_Deployments,'k-','LineWidth',4)
yscale('log')
xscale('log')
set(gca,'FontSize',20)
title('Mean FFT: Deployment #1 & #2 with Chunking')
grid on