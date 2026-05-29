clear all
close all
clc

script_dir = fileparts(mfilename('fullpath'));
data_folder = fullfile(script_dir,'matfiles');
load(fullfile(data_folder,'openloop_data_train.mat'),'y','u','t');

for n = 1:2
    if n==2
        load(fullfile(data_folder,'openloop_data_test.mat'),'y','u','t');
    end

    figure
    theme("light")
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
end

