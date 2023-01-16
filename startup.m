function startup()
format long

% plotting preferences
set(groot,'defaultAxesTickLabelInterpreter', 'latex')
set(groot,'defaultColorbarTickLabelInterpreter', 'latex')
set(groot,'defaultLegendInterpreter', 'latex')
set(groot,'defaultTextInterpreter', 'latex')
set(groot,'defaultAxesFontSize', 14)
set(groot,'defaultLegendFontSize', 14)

% adding current folder and subfolders to matlab path
addpath(genpath(pwd))
end