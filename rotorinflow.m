clear all;
close all;
clc;

speed = [10 20 40 50 60 80 100];
inflow_mean = zeros(1,2,length(speed));

for s = 1:length(speed)
    % Loading the Standalone TR data
    inflow_TR = 'C:\Users\hauti\heli-fw\MyCasesNew\t';
    %Loading the Standalone MR data from HeliUM
    inflow_MR = 'C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_';
    % Loading and concantenating along the 3rd dimensions the 2 rotor cases - nn & nb 
    inflow_2Rnn = 'C:\Users\hauti\heli-fwmrtr2\MyCases2rot\d';
    inflow_2Rnb = 'C:\Users\hauti\heli-fwmrtr2\MyCases2rot\d';
    filend = {'nn-nzt577', 'nb-nzt577'};
    count = 1;
    % Extracing the inflow values from the standalone TR and MR rotor cases.
    [~, inflownorm] = importINFLOW(speed(s), inflow_TR, 1,'');
    inflow_mean(count,:,s) = inflownorm;
    count = count + 1;
    % Extracing the inflow values from the standalone MR rotor case.
    % [inflowmat, inflownorm] = importINFLOW(speed, inflow_MR, 1,'');
    % inflow_mean(count) = inflownorm;
    % count = count + 1;
    % Extracing the inflow values from the two rotor case 'nn'.
    [~, inflownorm] = importINFLOW(speed(s), inflow_2Rnn, 2,'nn-nzt577');
    inflow_mean(count,:,s) = inflownorm;
    count = count + 1;
    % Extracing the inflow values from the two rotor case 'nb'.
    [~, inflownorm] = importINFLOW(speed(s), inflow_2Rnb, 2,'nb-nzt577');
    inflow_mean(count,:,s) = inflownorm;
end

%% Plotting the inflows
  
figure(1) % Plotting the inflow for the main rotor
inflowMR(1:2,:) = inflow_mean(2:3,1,:);
plot(speed, inflowMR(1,:), 'r', speed,inflowMR(2,:),'b')
xlabel('Velocity in kts')
ylabel('Inflow coefficient')
title('MAIN ROTOR')
legend('2 rotor - nn', '2 rotor - nb')

% figure(2) % Plotting the inflow for the tail rotor
inflowTR(1:3,:) = inflow_mean(1:3,2,:);
figure(2)
plot(speed,inflowTR(1,:),'r','LineWidth',1.6)
hold on;
plot(speed,inflowTR(2,:),'b','LineWidth',1.6)
plot(speed,inflowTR(3,:),'g','LineWidth',1.6)
xlabel('Velocity in kts')
ylabel('Inflow coefficient')
title('TAIL ROTOR')
legend('Standalone Rotor','2 rotor - nn', '2 rotor - nb')





