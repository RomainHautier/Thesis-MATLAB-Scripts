clear;
clc;
close all;
figcount = 1;

% This script takes in control angles and converts them to pedal inputs in
% inches and percentages.
saveDir = ['C:\Users\hauti\OneDrive - Imperial College London\Documents\London 23-24\Research Project\Result Figures\Inflow Plots\'];
%% Extracting the Free Wake Standalone data.

speed = 10:5:100;

for i = 1:length(speed)
    [timeflightmat, timeflight, t0rfw, t1cfw, t1sfw, b0r, b1c, b1s, timeflight_lastline] = importtimeflight(speed, 'C:\Users\hauti\heli-fwmrtr2\MyCasesNew-8turns\TR\t', '');
end

% Loading the TR collective results from Helium
t0rhel = struct2array(load("t0rHel.mat"));
exp_data = struct2array(load("exp_pedal_data.mat"));

%% Converting the Control Angles to Pedal Inputs

neutral_t0r_TR_deg = 15; % Neutral collective pitch angle of the TR (deg)
neutral_t0r_TR_in = 2.69; % Neutral collective pitch angle of the TR (inches)
pedal_range = 29.8; % deg

pedal_percentage = (pedal_range - t0rfw)/pedal_range*100;
pedal_percentage_hel = (pedal_range - t0rhel)/pedal_range*100;
pedal_percentage_tt = (pedal_range - 21.5711)/pedal_range*100

% t0r_TR_pedal = neutral_t0r_TR_deg + (t0rfw - neutral_t0r_TR_deg).*20.66/; % Tail rotor collective pitch input conversion from inches to degrees.
% t0r_TR = t0r_TR_pedal + flip(TTR.*S6 + BIASTR);


%% Creating a correction accounting for the change in inflow to calculate a new pedal collective.
om = 124.6;
rad = 1.676;
V = om*rad; % TR blade tip speed

inflowdat = struct2array(load("inflowTR2rot.mat"));
inflow2rot = inflowdat(2,:);
inflowhel = struct2array(load("inflowTRHel.mat"));

for i = 3:length(inflowhel)
    dT(i-2) = atand(inflow2rot(i-2)) - atand(inflowhel(i));
end

t0r2rot =  t0rhel(3:end) + dT';
pedal_percentage_2rot = (pedal_range - t0r2rot)/pedal_range*100;


%% Converstion trial
pedal_inches = 1.932667
pedal_degrees = 29.9 - pedal_inches*29.8/5.38
% t0r_TR_pedal = neutral_t0r_TR_deg + (t0r_helium_ouput - neutral_t0r_TR_in).*29.8/5.38; % Tail rotor collective pitch input conversion from inches to degrees.

%% Plotting Section
figure(figcount)
plot([1 5 speed],pedal_percentage_hel,'color',[0.9290 0.6940 0.1250],'LineWidth',1.6)
hold on;
scatter([1 5 speed],pedal_percentage_hel,20,[0.9290 0.6940 0.1250],'filled')
plot(speed,pedal_percentage,'color',[0 0.4470 0.7410],'LineWidth',1.6)
scatter(speed,pedal_percentage,20,[0 0.4470 0.7410],'filled')
plot(speed,pedal_percentage_2rot,'color',[0.6350 0.0780 0.1840],'LineWidth',1.6)
scatter(speed,pedal_percentage_2rot,20,[0.6350 0.0780 0.1840],'filled')
scatter(exp_data(:,1),exp_data(:,2),36,[0 0 0],'filled')
plot([0,100],[50 50],'color',[0 0 0],'LineWidth',1.6,'LineStyle','--')
xlabel('Forward Speed $V$ (kts)','Interpreter','latex','FontSize',14)
ylabel('Pedal Position $\delta_{ped}$ ($\%$)','Interpreter','latex','FontSize',14)
legend('HeliUM','','Standalone Free Wake','','Full Inflow Correction','','Experimental Data','Pedal Neutrality','Interpreter','latex','Fontsize',12)
xticks(0:10:100)
yticks(0:10:100)
allAxes = findall(gcf,'Type','axes'); % Find all axes in the figure
set(allAxes, 'XLim', [10 100]);         % Set common xlim
set(allAxes, 'XGrid', 'on', 'YGrid', 'on');
set(allAxes, 'XTick', 10:10:100) %sets the x-tick marks at intervals of 2 from 0 to 10 for all tiles.
set(allAxes, 'XMinorTick', 'on')
set(allAxes, 'YMinorTick', 'on') %enable minor ticks for both x and y axes.

ylim([0 100])


%% Plotting the Tail Rotor Collective vs V for all 3 models used.

