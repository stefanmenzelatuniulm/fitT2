close all;
clear all;
clc;

%-------------SETTINGS-------------

echoSpacing = 25.8; %ms;
numberOfEchoes = 150;
path = "C:\Users\menze\Desktop\Matlab\MR_Data\2024_03_22_1";
chemicalSpecies = "Lactate"; %Name(s) of the chemical species
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot

%-------------END OF SETTINGS-------------

load(path+"\data.mat");
load(path+"\par.mat");

%Data_raw = importdata(path+"\sp_21032024_1520298_98_2_wip_t1w_tseV4.raw");
dim = size(Data);

%Take all phase encoding lines and concatenate them horizontally
Data = reshape(Data, [1 dim(2)*dim(1) dim(3) dim(4)]);

%FFT for all coil elements separatly
for k = 1:dim(4)
    Data(:, 1, 1, k) = fft(Data(:, 1, 1, k));
end

%Sum over coil elements
Data = abs(Data);
Data = squeeze(sum(Data, 4));
Data = double(Data/max(Data));

fig = figure('WindowState', 'maximized');
X = transpose(linspace(1, echoSpacing*numberOfEchoes, length(Data)));
plot(X, Data, "+");
fitfunction="C+M0*exp(-x/T2)";
coeffs=["C" "M0" "T2"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-1 0 0.1], 'Upper', [1 1 inf], 'StartPoint', [0 0 1000]);
fttype = fittype(fitfunction, coefficients=coeffs);
ft=fit(X, transpose(Data), fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_2", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(X)]);
legend("Signal", "Fit with $$C+M_0 e^{-\frac{t}{T_2}}$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("T2 decay during TSE of "+chemicalSpecies);
xlabel("$t$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_TSEDecay.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay.svg");

close all;