clear all;
close all;
clc;

%-------------SETTINGS-------------

%tau = exp(linspace(log(10),log(10000),12)); %ms;
%tau = [0 1 10 100 250 527.3 1053.5 1579.8 2106.1 2632.3 3158.6 3684.8 4211.1 4737.4 5263.6 5789.9 6316.2 6842.4 7368.7 7894.9 8421.2 8947.5 9473.7 10000];
%tau = [0 1 10 100 250 527.3 1053.5 1579.8];
%tau = [0 2632.3 3158.6 100 250 527.3 1053.5 1579.8 2106.1 4211.1 5263.6 6316.2 7368.7 8421.2 9473.7 10000];
tau = [0 100 250 527.3 1053.5 1579.8 2106.1 2632.3 3158.6 4211.1 5263.6 6316.2 7368.7 8421.2 9473.7 10000];
%tau = [0 100 250 527.3 1579.8 2106.1 2632.3 3158.6 4211.1 5263.6 6316.2 7368.7 8421.2 9473.7 10000];
folderName = "C:\Users\Stefan Menzel\Desktop\Matlab\MR_Data\2024_04_19\Alanin\T1Ala";
chemicalSpecies = "Alanine"; %Name(s) of the chemical species
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot
outlierIndices = [5];

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    dim = size(Data);
    if mod(dim(2), 2) == 0
        Data1 = Data(k, ceil(end/2), :);
        Data2 = Data(k, ceil(end/2)+1, :);
        Data(k, 1, :) = (Data1+Data2)/2;
    else
        Data(k, 1, :) = Data(k, ceil(end/2), :);
    end
    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 

% [m, minIndex] = min(Data); %Remove bounce plot
% for k = 1:minIndex-1
%     Data(k) = -Data(k);
% end

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-1 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center_flipped_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center_flipped_absLast.svg");

close all;

%Alternatively: Take maximum of signal

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max_flipped_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max_flipped_absLast.svg");

close all;

%Alternatively: Take mean in readout direction and max in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max_flipped_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max_flipped_absLast.svg");

close all;

%Alternatively: Take mean in readout direction and center in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
dim = size(Data);
clear Data2;

for k = 1:length(dFolders)

    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center_flipped_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center_flipped_absLast.svg");

close all;

%Alternatively: Take mean in readout direction and mean in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
Data = mean(Data, 3);
dim = size(Data);
clear Data2;

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean_flipped_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean_flipped_absLast.svg");

close all;

%Alternatively: Assume that odd and even averages are not flipped

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data = sum(Data2, 13);
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    dim = size(Data);
    if mod(dim(2), 2) == 0
        Data1 = Data(k, ceil(end/2), :);
        Data2 = Data(k, ceil(end/2)+1, :);
        Data(k, 1, :) = (Data1+Data2)/2;
    else
        Data(k, 1, :) = Data(k, ceil(end/2), :);
    end
    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 

% [m, minIndex] = min(Data); %Remove bounce plot
% for k = 1:minIndex-1
%     Data(k) = -Data(k);
% end

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-1 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center_absLast.svg");

close all;

%Alternatively: Take maximum of signal

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data = sum(Data2, 13);
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max_absLast.svg");

close all;

%Alternatively: Take mean in readout direction and max in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data = sum(Data2, 13);
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max_absLast.svg");

close all;

%Alternatively: Take mean in readout direction and center in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data = sum(Data2, 13);
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
dim = size(Data);
clear Data2;

for k = 1:length(dFolders)

    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center_absLast.svg");

close all;

%Alternatively: Take mean in readout direction and mean in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data = sum(Data2, 13);
Data = double(abs(Data)); %FID
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
Data = mean(Data, 3);
dim = size(Data);
clear Data2;

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean_absLast.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean_absLast.svg");

close all;

%Alternatively: abs first

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    dim = size(Data);
    if mod(dim(2), 2) == 0
        Data1 = Data(k, ceil(end/2), :);
        Data2 = Data(k, ceil(end/2)+1, :);
        Data(k, 1, :) = (Data1+Data2)/2;
    else
        Data(k, 1, :) = Data(k, ceil(end/2), :);
    end
    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 

% [m, minIndex] = min(Data); %Remove bounce plot
% for k = 1:minIndex-1
%     Data(k) = -Data(k);
% end

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-1 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center_flipped.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center_flipped.svg");

close all;

%Alternatively: Take maximum of signal

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max_flipped.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max_flipped.svg");

close all;

%Alternatively: Take mean in readout direction and max in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max_flipped.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max_flipped.svg");

close all;

%Alternatively: Take mean in readout direction and center in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
dim = size(Data);
clear Data2;

for k = 1:length(dFolders)

    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center_flipped.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center_flipped.svg");

close all;

%Alternatively: Take mean in readout direction and mean in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data_odd = sum(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 1:2:end), 13);
Data_even = sum(flip(Data2(:, :, :, :, :, :, :, :, :, :, :, :, 2:2:end), 2), 13);
Data = Data_odd+Data_even;
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
Data = mean(Data, 3);
dim = size(Data);
clear Data2;

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean_flipped.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean_flipped.svg");

close all;

%Alternatively: Assume that odd and even averages are not flipped

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data = sum(Data2, 13);
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    dim = size(Data);
    if mod(dim(2), 2) == 0
        Data1 = Data(k, ceil(end/2), :);
        Data2 = Data(k, ceil(end/2)+1, :);
        Data(k, 1, :) = (Data1+Data2)/2;
    else
        Data(k, 1, :) = Data(k, ceil(end/2), :);
    end
    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 

% [m, minIndex] = min(Data); %Remove bounce plot
% for k = 1:minIndex-1
%     Data(k) = -Data(k);
% end

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-1 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_center.svg");

close all;

%Alternatively: Take maximum of signal

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data = sum(Data2, 13);
Data = squeeze(sum(Data, 5)); %Sum over coils
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_max.svg");

close all;

%Alternatively: Take mean in readout direction and max in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data = sum(Data2, 13);
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
clear Data2;

for k = 1:length(dFolders)

    Data(k, 1, 1) = max(Data(k, :, :), [], "all");

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_max.svg");

close all;

%Alternatively: Take mean in readout direction and center in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data = sum(Data2, 13);
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
dim = size(Data);
clear Data2;

for k = 1:length(dFolders)

    if mod(dim(3), 2) == 0
        Data1 = Data(k, :, ceil(end/2));
        Data2 = Data(k, :, ceil(end/2)+1);
        Data(k, :, 1) = (Data1+Data2)/2;
    else
        Data(k, :, 1) = Data(k, :, ceil(end/2));
    end

end

Data = Data(:, 1, 1); 
Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_center.svg");

close all;

%Alternatively: Take mean in readout direction and mean in phase encoding
%direction

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
    Data2(k, :, :, :, :, :, :, :, :, :, :, :, :) = Data;
end
Data2 = double(abs(Data2)); %FID
Data = sum(Data2, 13);
Data = squeeze(sum(Data, 5)); %Sum over coils
Data = mean(Data, 2);
Data = mean(Data, 3);
dim = size(Data);
clear Data2;

Data = Data/max(Data, [], "all"); %Normalize

fig = figure('WindowState', 'maximized');
plot(tau, Data, "+");
fitfunction="abs(M0*(1-C*exp(-x/T1)))";
coeffs=["M0" "C" "T1"];
options=fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0 0.1], 'Upper', [1 inf inf], 'StartPoint', [1 2 9000]);
fttype = fittype(fitfunction, coefficients=coeffs);
fitData = Data(setdiff(1:end, outlierIndices, "sorted"));
fitTau = tau(setdiff(1:end, outlierIndices, "sorted"));
ft=fit(transpose(fitTau), fitData, fttype, options);
coeffvals = coeffvalues(ft);
ci = confint(ft, 0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)', "T_1", coeffvals(3), ci(:, 3));
annotation('textbox', [0.53+annotationXOffset 0.69 0.2 0.2], 'String', ['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")], 'EdgeColor', 'none', "FitBoxToText", "on", "Color", "r", "FontSize", 8);
hold on;
ax = gca;
plot(ax, ft, "r");
xlim([0 max(tau)]);
legend("Signal", "Fit with $$|M_0(1-C e^{-\frac{\tau}{T_1}})|$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("$T_1$ Relaxation of "+chemicalSpecies+" during IR", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlabel("$\tau$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean.fig");
saveas(fig, folderName+"\"+num2str(subFolders(1,1))+chemicalSpecies+"_IRDecay_mean_mean.svg");

close all;


