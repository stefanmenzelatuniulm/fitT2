close all;
clear all;
clc;

%-------------SETTINGS-------------

echoSpacing = 25.8; %ms;
path = "C:\Users\Stefan Menzel\Desktop\Matlab\MR_Data\2024_04_19\T2Ala\197";
chemicalSpecies = "Alanine"; %Name(s) of the chemical species
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot

%-------------END OF SETTINGS-------------

%Better to take abs() first before summing over averages? -> Closer to
%Julians result. Moreover, fitting model assumes that the signal decays towards
%a constant noise floor instead of 0

%Reco-PC returns complex signal (not spectrum)
load(path+"\data.mat");
load(path+"\par.mat");

dim = size(Data);

%Data = sum(Data, 12); 
%Sum over averages
%CARE: ODD AND EVEN AVERAGES ARE FLIPPED
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 4)); %Sum over coils

%Echo tops are centers of k-space lines in readout direction
if mod(dim(1), 2) == 0
    Data1 = Data(ceil(end/2), :, :, :);
    Data2 = Data(ceil(end/2)+1, :, :, :);
    Data = (Data1+Data2)/2;
else
    Data = Data(ceil(end/2), :, :, :);
end

Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_absLast.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_absLast.svg");

%ALTERNATIVELY: Echo Tops have approx. all the signal, rest of readout
%dimension has approx. no signal if T2* << T2 -> mean over readout gives
%signal at echo tops

load(path+"\data.mat");
load(path+"\par.mat");

%Data = sum(Data, 12); 
%Sum over averages
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 4)); %Sum over coils
Data = mean(Data, 1); %Mean over readout direction
Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_mean_absLast.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_mean_absLast.svg");

close all;

%ALTERNATIVELY: Echo Tops are located at the max signal for every row
load(path+"\data.mat");
load(path+"\par.mat");

%Data = sum(Data, 12); 
%Sum over averages
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 4)); %Sum over coils
Data = max(Data, [], 1); 
Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_max_absLast.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_max_absLast.svg");

close all;

%ALTERNATIVELY: Echo Tops are located at the row with the overall max signal
load(path+"\data.mat");
load(path+"\par.mat");

%Data = sum(Data, 12); 
%Sum over averages
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 4)); %Sum over coils
m = max(Data, [], "all"); 
[r, c] = find(abs(Data-m)<10^(-6));
Data = Data(r(1), :);
Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_max2_absLast.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_max2_absLast.svg");

close all;

%ALTERNATIVELY: abs() first
%Reco-PC returns complex signal (not spectrum)
load(path+"\data.mat");
load(path+"\par.mat");

dim = size(Data);

%Data = sum(Data, 12); 
%Sum over averages
%CARE: ODD AND EVEN AVERAGES ARE FLIPPED
Data = double(abs(Data)); %FID
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 4)); %Sum over coils

%Echo tops are centers of k-space lines in readout direction
if mod(dim(1), 2) == 0
    Data1 = Data(ceil(end/2), :, :, :);
    Data2 = Data(ceil(end/2)+1, :, :, :);
    Data = (Data1+Data2)/2;
else
    Data = Data(ceil(end/2), :, :, :);
end

Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

%ALTERNATIVELY: Echo Tops have approx. all the signal, rest of readout
%dimension has approx. no signal if T2* << T2 -> mean over readout gives
%signal at echo tops

load(path+"\data.mat");
load(path+"\par.mat");

%Data = sum(Data, 12); 
%Sum over averages
Data = double(abs(Data)); %FID
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 4)); %Sum over coils
Data = mean(Data, 1); %Mean over readout direction
Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_mean.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_mean.svg");

close all;

%ALTERNATIVELY: Echo Tops are located at the max signal for every row
load(path+"\data.mat");
load(path+"\par.mat");

%Data = sum(Data, 12); 
%Sum over averages
Data = double(abs(Data)); %FID
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 4)); %Sum over coils
Data = max(Data, [], 1); 
Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_max.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_max.svg");

close all;

%ALTERNATIVELY: Echo Tops are located at the row with the overall max signal
load(path+"\data.mat");
load(path+"\par.mat");

%Data = sum(Data, 12); 
%Sum over averages
Data = double(abs(Data)); %FID
Data_odd = sum(Data(:, :, :, :, :, :, :, :, :, :, :, 1:2:end), 12);
Data_even = sum(flip(Data(:, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 12);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 4)); %Sum over coils
m = max(Data, [], "all"); 
[r, c] = find(abs(Data-m)<10^(-6));
Data = Data(r(1), :);
Data = transpose(Data/max(Data, [], "all")); %Normalize

%Determine echo placements. Each echo corresponds to 1 line in k-space
dim = size(Data);
numberOfPhaseEncodingSteps = dim(1);
X = transpose(echoSpacing:echoSpacing:echoSpacing*numberOfPhaseEncodingSteps);

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

saveas(fig, path+chemicalSpecies+"_TSEDecay_max2.fig");
saveas(fig, path+chemicalSpecies+"_TSEDecay_max2.svg");

close all;