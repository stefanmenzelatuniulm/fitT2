close all;
clear all;
clc;

%-------------SETTINGS-------------

echoSpacing = 19.7; %ms;
path = "C:\Users\Stefan Menzel\Desktop\Matlab\MR_Data\2024_03_26\T2Bic2\157";
chemicalSpecies = "Bicarbonate"; %Name(s) of the chemical species
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot

%-------------END OF SETTINGS-------------

%Reco-PC returns complex signal (not spectrum)
load(path+"\data.mat");
load(path+"\par.mat");

Data = abs(Data); %FID
Data = mean(Data, 1); %Mean over readout direction / for every k-space-line
Data = permute(Data, [2 4 1 3]); %dimensions k-space-lines, coils
Data = sum(Data, 2); %Sum over coils
Data = double(Data/max(Data)); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing/2:echoSpacing:echoSpacing*(2*numberOfPhaseEncodingSteps-1)/2);

fig = figure('WindowState', 'maximized');
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