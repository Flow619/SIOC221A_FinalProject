%% Trenton_SIOC221A_Drifter_Analysis_Chunking
%
%% 
close all
clear all
clc

%%
addpath('C:\Users\Trenton\Documents\GitHub\SIOC221A\HW4')
addpath('C:\Users\Trenton\Documents\GitHub\SIOC221A\HW3')

%% User Defined
Deployment_Num = 1; % Either #1 or #2
dt = 0.2; % sec;

Save_Flag = "OFF";

%% Load Data
if Deployment_Num == 1
    FilePath = 'C:\Users\Trenton\Documents\GitHub\SIOC221A_FinalProject\drifter\data_dep1';
elseif Deployment_Num == 2
    FilePath = 'C:\Users\Trenton\Documents\GitHub\SIOC221A_FinalProject\drifter\data_dep2';
end

cd(FilePath)
CSV_Files = dir;
CSV_Files = {CSV_Files(3:end-1).name}'; % pull list of .csv names

Hanning = hann(600); % hanning window

if Deployment_Num == 1
    for i = 3:length(CSV_Files)-2 % ***Only pull .csv files when drifter was in the water***

        Data = readtable(CSV_Files{i});
        millis{i-2} = Data.millis;

        % millg to g
        AccZ{i-2} = Data.accZ/1000;
        AccY{i-2} = Data.accY/1000;
        AccX{i-2} = Data.accX/1000;

        % demean
        AccZ{i-2} = AccZ{i-2} - mean(AccZ{i-2});
        AccY{i-2} = AccY{i-2} - mean(AccY{i-2});
        AccX{i-2} = AccX{i-2} - mean(AccX{i-2});

        % FFT (with hanning and demeaning)
        [P_accZ_Chunk1{i-2},freq] = MySpectrum( (Hanning.*AccZ{i-2}(1:600))',dt,"OFF");
        [P_accY_Chunk1{i-2},freq] = MySpectrum( (Hanning.*AccY{i-2}(1:600))',dt,"OFF");
        [P_accX_Chunk1{i-2},freq] = MySpectrum( (Hanning.*AccX{i-2}(1:600))',dt,"OFF");

        [P_accZ_Chunk2{i-2},freq] = MySpectrum( (Hanning.*AccZ{i-2}(151:750))',dt,"OFF");
        [P_accY_Chunk2{i-2},freq] = MySpectrum( (Hanning.*AccY{i-2}(151:750))',dt,"OFF");
        [P_accX_Chunk2{i-2},freq] = MySpectrum( (Hanning.*AccX{i-2}(151:750))',dt,"OFF");
        
        [P_accZ_Chunk3{i-2},freq] = MySpectrum( (Hanning.*AccZ{i-2}(301:900))',dt,"OFF");
        [P_accY_Chunk3{i-2},freq] = MySpectrum( (Hanning.*AccY{i-2}(301:900))',dt,"OFF");
        [P_accX_Chunk3{i-2},freq] = MySpectrum( (Hanning.*AccX{i-2}(301:900))',dt,"OFF");

        % Normalize for Hanning
        P_accZ_Chunk1{i-2} = P_accZ_Chunk1{i-2}/ mean(Hanning.^2);
        P_accZ_Chunk2{i-2} = P_accZ_Chunk2{i-2}/ mean(Hanning.^2);                        
        P_accZ_Chunk3{i-2} = P_accZ_Chunk3{i-2}/ mean(Hanning.^2);

        % Can repeat this process for X,Y,Z if desired.
        % P_accY{i-2} = P_accY{i-2}/ mean(Hanning.^2);
        % P_accX{i-2} = P_accX{i-2}/ mean(Hanning.^2);

    end
elseif Deployment_Num == 2
    for i = 13:length(CSV_Files)-1 % ***Only pull .csv files when drifter was in the water***

        Data = readtable(CSV_Files{i});
        millis{i-12} = Data.millis;

        % millg to g
        AccZ{i-12} = Data.accZ/1000;
        AccY{i-12} = Data.accY/1000;
        AccX{i-12} = Data.accX/1000;

        % Demean
        AccZ{i-12} = AccZ{i-12} - mean(AccZ{i-12});
        AccY{i-12} = AccY{i-12} - mean(AccY{i-12});
        AccX{i-12} = AccX{i-12} - mean(AccX{i-12});

        % FFT (with hanning and demeaning)
        [P_accZ_Chunk1{i-12},freq] = MySpectrum((Hanning.*AccZ{i-12}(1:600))',dt,"OFF");
        [P_accY_Chunk1{i-12},freq] = MySpectrum((Hanning.*AccY{i-12}(1:600))',dt,"OFF");
        [P_accX_Chunk1{i-12},freq] = MySpectrum((Hanning.*AccX{i-12}(1:600))',dt,"OFF");

        [P_accZ_Chunk2{i-12},freq] = MySpectrum((Hanning.*AccZ{i-12}(151:750))',dt,"OFF");
        [P_accY_Chunk2{i-12},freq] = MySpectrum((Hanning.*AccY{i-12}(151:750))',dt,"OFF");
        [P_accX_Chunk2{i-12},freq] = MySpectrum((Hanning.*AccX{i-12}(151:750))',dt,"OFF");

        [P_accZ_Chunk3{i-12},freq] = MySpectrum((Hanning.*AccZ{i-12}(301:900))',dt,"OFF");
        [P_accY_Chunk3{i-12},freq] = MySpectrum((Hanning.*AccY{i-12}(301:900))',dt,"OFF");
        [P_accX_Chunk3{i-12},freq] = MySpectrum((Hanning.*AccX{i-12}(301:900))',dt,"OFF");


        % Normalize for Hanning
        P_accZ_Chunk1{i-12} = P_accZ_Chunk1{i-12}/ mean(Hanning.^2);
        P_accZ_Chunk2{i-12} = P_accZ_Chunk2{i-12}/ mean(Hanning.^2);
        P_accZ_Chunk3{i-12} = P_accZ_Chunk3{i-12}/ mean(Hanning.^2);
        
        % Can repeat this process for X,Y,Z if desired.
        % P_accY{i-12} = P_accY{i-12}/ mean(Hanning.^2);
        % P_accX{i-12}= P_accX{i-12}/ mean(Hanning.^2);
    end
end

%% Calculate Mean FFT
Mean_P = mean([P_accZ_Chunk1{:},P_accZ_Chunk2{:},P_accZ_Chunk3{:}],2);

%% Plot FFTs

figure('units','normalized','outerposition',[0 0 1 1])
hold on

for ii = 1:size(P_accZ_Chunk1,2)
    plot(freq',P_accZ_Chunk1{ii},'linewidth',1) %Plot each FFT
    plot(freq',P_accZ_Chunk2{ii},'linewidth',1) %Plot each FFT
    plot(freq',P_accZ_Chunk3{ii},'linewidth',1) %Plot each FFT

end

yscale('log')
xscale('log')
ylim([10^-9 10^0])

plot(freq,Mean_P,'k.-','LineWidth',4,'HandleVisibility','Off')
set(gca,'FontSize',18)

title(['Deployment #',num2str(Deployment_Num), ' with Chunking' ],'FontSize',20)

%%  Save FFT Data
if strcmp(Save_Flag,"ON")
    cd('C:\Users\Trenton\Documents\GitHub\SIOC221A_FinalProject\Trenton_Data_Code\FIGURES')
    if Deployment_Num == 1
        Deploy1 = [P_accZ_Chunk1{:},P_accZ_Chunk2{:},P_accZ_Chunk3{:}];
        save('FFT_Deploy1_Chunking','Deploy1',"freq")
        exportgraphics(gcf,'FFT_Deployment1_Chunking.jpg','Resolution',300)
    elseif Deployment_Num == 2
        Deploy2 = [P_accZ_Chunk1{:},P_accZ_Chunk2{:},P_accZ_Chunk3{:}];
        save('FFT_Deploy2_Chunking','Deploy2')
        exportgraphics(gcf,'FFT_Deployment2_Chunking.jpg','Resolution',300)
    end
end
