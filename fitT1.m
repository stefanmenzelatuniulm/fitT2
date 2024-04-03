clear all;
close all;
clc;

%-------------SETTINGS-------------

tau = exp(linspace(log(10),log(10000),12)); %ms;
folderName = "C:\Users\Stefan Menzel\Desktop\Matlab\MR_Data\2024_03_28\T1Lac";
chemicalSpecies = "Lactate"; %Name(s) of the chemical species
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot

%-------------END OF SETTINGS-------------

d = dir(folderName);
dFolders = d([d(:).isdir]);
dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));
subFolders = zeros(1, length(dFolders));
Data = [];
for k = 1:length(dFolders)
    subFolders(1, k) = string(dFolders(k).name);
end
load(folderName+"\"+subFolders(1, k)+"\data.mat");
dim = size(Data);
Data2 = zeros([length(dFolders), dim]);
for k = 1:length(dFolders)
    load(folderName+"\"+subFolders(1, k)+"\data.mat");
    Data2(k, :, :, :, :) = Data;
end
Data = abs(Data2);
clear Data2;

for k = 1:length(dFolders)

    dim = size(Data(k, :, :, :, :));
    if mod(dim(2), 2) == 0
        Data1 = Data(k, ceil(end/2), :, :, :);
        Data2 = Data(k, ceil(end/2)+1, :, :, :);
        Data(k, 1, :, :, :) = (Data1+Data2)/2;
    else
        Data(k, 1, :, :, :) = Data(k, ceil(end/2), :, :, :);
    end
    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2), :, :);
        Data2 = Data(k, :, ceil(end/2)+1, :, :);
        Data(k, :, 1, :, :) = (Data1+Data2)/2;
    else
        Data(k, :, 1, :, :) = Data(k, :, ceil(end/2), :, :);
    end

end

Data = permute(Data(:, 1, 1, :, :), [1 5 2 3 4]); %indices Tau, Coil
Data = sum(Data, 2); %Sum over coils

% [m, minIndex] = min(Data); %Remove bounce plot
% for k = 1:minIndex-1
%     Data(k) = -Data(k);
% end

Data = double(Data/max(Data)); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 2 inf], 'StartPoint', [1 2 1000]);
fttype = fittype(fitfunction, coefficients=coeffs);
ft=fit(transpose(tau), Data, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("T1 relaxation of "+chemicalSpecies+" during IR");
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+chemicalSpecies+"_IRDecay.fig");
saveas(fig, folderName+chemicalSpecies+"_IRDecay.svg");

close all;

%Alternatively: Sum over whole signal: Every signal point ~ to T1
%relaxation -> therefore sum over whole signal ~ to T1 relaxation

d = dir(folderName);
dFolders = d([d(:).isdir]);
dFolders = dFolders(~ismember({dFolders(:).name},{'.','..'}));
subFolders = zeros(1, length(dFolders));
Data = [];
for k = 1:length(dFolders)
    subFolders(1, k) = string(dFolders(k).name);
end
load(folderName+"\"+subFolders(1, k)+"\data.mat");
dim = size(Data);
Data2 = zeros([length(dFolders), dim]);
for k = 1:length(dFolders)
    load(folderName+"\"+subFolders(1, k)+"\data.mat");
    Data2(k, :, :, :, :) = Data;
end
Data = abs(Data2);
clear Data2;

for k = 1:length(dFolders)

    dim = size(Data(k, :, :, :, :));
    Data(k, 1, :, :, :) = sum(Data(k, :, :, :, :), 2);
    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2), :, :);
        Data2 = Data(k, :, ceil(end/2)+1, :, :);
        Data(k, :, 1, :, :) = (Data1+Data2)/2;
    else
        Data(k, :, 1, :, :) = Data(k, :, ceil(end/2), :, :);
    end

end

Data = abs(Data(:, 1, 1, :, :));
Data = permute(Data(:, 1, 1, :, :), [1 5 2 3 4]); %indices Tau, Coil
Data = sum(Data, 2); %Sum over coils
% [m, minIndex] = min(Data); %Remove bounce plot
% for k = 1:minIndex-1
%     Data(k) = -Data(k);
% end
Data = double(Data/max(Data)); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 2 inf], 'StartPoint', [1 2 1000]);
fttype = fittype(fitfunction, coefficients=coeffs);
ft=fit(transpose(tau), Data, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("T1 relaxation of "+chemicalSpecies+" during IR");
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+chemicalSpecies+"_IRDecay_mean.fig");
saveas(fig, folderName+chemicalSpecies+"_IRDecay_mean.svg");

close all;