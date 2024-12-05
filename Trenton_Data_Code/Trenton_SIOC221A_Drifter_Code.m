%SIOC221A Final Project Data Analysis Code
%
% Data is kind of fucked. treat each .csv file indep and then apply hanning
% window. Take the average.
%% 
close all;
clear all;
clc;

addpath('C:\Users\Trenton\Documents\GitHub\SIOC221A\HW4')
addpath('C:\Users\Trenton\Documents\GitHub\SIOC221A\HW3')

%% User Defined
Deployment_Num = 1; %Either #1 or #2

%% Load Data
if Deployment_Num == 1
    FilePath = 'C:\Users\Trenton\Documents\GitHub\SIOC221A_FinalProject\drifter\data_dep1'
elseif Deployment_Num == 2
    FilePath = 'C:\Users\Trenton\Documents\GitHub\SIOC221A_FinalProject\drifter\data_dep2'
end

cd(FilePath)
CSV_Files = dir;
CSV_Files = {CSV_Files(3:end-1).name}';

if Deployment_Num == 1
    for i = 3:length(CSV_Files)-1
        if i == 3
            Data = readtable(CSV_Files{i});
        else
            Temp = readtable(CSV_Files{i});
            Data = [Data;Temp];
        end
    end
elseif Deployment_Num == 2
    for i = 12:length(CSV_Files)
        if i == 12
            Data = readtable(CSV_Files{i});
        else
            Temp = readtable(CSV_Files{i});
            Data = [Data;Temp];
        end
    end
end

%% Truncate usable timeseries
if Deployment_Num == 1

x = Data.accZ(1750:10000)';
dt = 200/1000;
Figure_Flag = "ON"
%% Fourier Analysis Full Time Series
[P,freq] = MySpectrum(x,dt,Figure_Flag)

Num_windows = 10;
[Mean_P,freq,Error_Bars] = SpectralFunction_2(x,dt,Num_windows)

hold on
plot(freq,Mean_P, 'r','linewidth',3)

end
%%