close all;
clear all;
clc;

%-------------SETTINGS-------------

echoSpacing = 25.8; %ms;
numberOfPhaseEncodingSteps = 155; %number of echo spacings is numberOfPhaseEncodingSteps-1
path = "C:\Users\Stefan Menzel\Desktop\Matlab\MR_Data\2024_03_22\T2Lac\114";
chemicalSpecies = "Lactate"; %Name(s) of the chemical species
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot

%-------------END OF SETTINGS-------------

%Reco-PC returns complex signal (not spectrum)
load(path+"\data.mat");
load(path+"\par.mat");

dim = size(Data);

%Take all phase encoding lines and concatenate them horizontally
Data = reshape(Data, [1 dim(2)*dim(1) dim(3) dim(4)]); 
Data = permute(Data, [2 4 1 3]); %dimensions samples, coils

%For every sample, sum absolute values from all coils and normalize
Data = abs(Data);
Data = sum(Data, 2);
Data = double(Data/max(Data));

fig = figure('WindowState', 'maximized');
X = transpose(linspace(1, echoSpacing*(numberOfPhaseEncodingSteps-1), length(Data)));
plot(X, Data, "+");
fitfunction="C+M0*exp(-x/T2)";
coeffs=["C" "M0" "T2"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-1 0 0.1], 'Upper', [1 1 inf], 'StartPoint', [0 0 1000]);
fttype = fittype(fitfunction, coefficients=coeffs);
ft=fit(X, Data, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_2", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(X)]);
legend("Signal", "Fit with $$C+M_0 e^{-\frac{t}{T_2}}$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("T2 decay of "+chemicalSpecies+" during TSE");
xlabel("$t$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_TSEDecay.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay.svg");

close all;