clear all
close all
clc

script_dir = fileparts(mfilename('fullpath'));
data_folder = fullfile(script_dir,'matfiles');
load(fullfile(data_folder,'openloop_data_1.mat'),'y','u','t');
load(fullfile(data_folder,'openloop_data_2.mat'),'y','u','t');
load(fullfile(data_folder,'openloop_data_3.mat'),'y','u','t');

figure
theme(figure,"light")
set(gcf, 'DefaultAxesFontSize', 24) 
set(gcf, 'DefaultTextFontSize', 24)
subplot(2,1,1), hold on, grid on   
plot(t,y(1,:),'.','MarkerSize',12)
plot(t,y(2,:),'.','MarkerSize',12)
legend('Temperature 1','Temperature 2','Location','best')
xlabel('Time [s]')
ylabel('Temperature [°C]')
subplot(2,1,2), hold on, grid on   
width = 5;
stairs(t,u(1,:),'LineWidth',width)
stairs(t,u(2,:),'LineWidth',width)
legend('Heater 1','Heater 2','Location','best')
xlabel('Time [s]')
ylabel('Heater control [%]')
ylim([0 100]);

%%
load(fullfile(data_folder,'model','MSE.mat'))

% Create bar graph
figure;
theme(figure,"light")
set(gcf, 'DefaultAxesFontSize', 24) 
set(gcf, 'DefaultTextFontSize', 24)
bar(MSE, 'FaceColor', [0.529, 0.808, 0.980]); % Sky blue color
xlabel('Order n');
ylabel('MSE');
title('MSE by Order n');
grid on;