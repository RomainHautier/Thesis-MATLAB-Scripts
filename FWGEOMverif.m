clear;
close all;
clc;

%% FWGEOM.DAT concatenation validity verification for 2 rotor cases.
% This function is dedicated to the verification of the validity of the
% initial wake geometry fed to the mfw solver for a 2 rotor case. To do so,
% the geometry of both MR & TR wakes are plotted for a single speed and
% compared to the original files they were pulled from. 

% Setting up the constants required for subsequent plotting
turns = 8;
dz = deg2rad(10);
range_of_speed = 10:5:100;
range_of_speed = [10 25 40 55 70 100];
rad = [8.18 1.676];
nzt = turns*2*pi/dz+1;
figcount = 1;

for i = 1:length(range_of_speed)
    
    % Opening a new figure & Setting up the tiled layout.
    figure(figcount)
    t = tiledlayout(1,2);
    nexttile
    
    % Plotting the wake geometry saved in the separate MR & TR freewake runs.
    [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(range_of_speed(i)) '\FWGEOM.DAT'],turns);
    % [pcx, pcy, pcz] = wakefilamentplot(1, ['C:\Users\hauti\heli-fwmrtr2\MyCasesNew-' num2str(turns) 'turns\MR\t' num2str(range_of_speed(i)) '\FWGEOM.DAT'],turns);
    plot3(pcx(36,1:nzt,1)*rad(1),pcy(36,1:nzt,1)*rad(1),pcz(36,1:nzt,1)*rad(1),'LineWidth',1.5)
    hold on
    % [pcx, pcy, pcz] = wakefilamentplot(1, ['C:\Users\hauti\heli-fwmrtr2\MyCasesNew-' num2str(turns) 'turns\TR\t' num2str(range_of_speed(i)) '\FWGEOM.DAT'],turns);
    % plot3(pcx(36,1:nzt,1)*rad(2),pcy(36,1:nzt,1)*rad(2),pcz(36,1:nzt,1)*rad(2),'LineWidth',1.5)
    plot3(pcx(36,1:nzt,2)*rad(2),pcy(36,1:nzt,2)*rad(2),pcz(36,1:nzt,2)*rad(2),'LineWidth',1.5)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Separate MR & TR Run Solutions')
    % view(0,0)
    nexttile
    
    % Plotting the wake geometry saved in the concatenated 2 rotor case.

    [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(range_of_speed(i)) '\IWGEOM.DATA'],turns);
    plot3(pcx(36,1:nzt,1)*rad(1),pcy(36,1:nzt,1)*rad(1),pcz(36,1:nzt,1)*rad(1),'LineWidth',1.5)
    hold on
    plot3(pcx(36,1:nzt,2)*rad(2),pcy(36,1:nzt,2)*rad(2),pcz(36,1:nzt,2)*rad(2),'LineWidth',1.5)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Concatenated MR & TR Run Solutions - IWGEOM.DATA')
    
    % Tiled layout attributes.
    title(t,['Freewake - 2 Rotor Case Starting Wake Solution at ' num2str(range_of_speed(i)) ' kts - ' num2str(turns) ' turns'])
    legend('Main Rotor', 'Tail Rotor')
    % view(0,0)
    figcount = figcount + 1;
end


% %% Checking the Similiarities in Solution between the MFW and MRTR scripts for a single TR rotor.
% figure(figcount)
% [pcx, pcy, pcz] = wakefilamentplot(1, ['C:\Users\hauti\heli-fw\MyCasesNew-' num2str(8) 'turns\TR\t' num2str(10) 'initial-mfw\FWGEOM.DAT'],8);
% plot3(pcx(36,1:289,1)*rad(2),pcy(36,1:289,1)*rad(2),pcz(36,1:289,1)*rad(2),'LineWidth',1.5)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('10kts - Case Ran with MFW Code')
% figcount = figcount + 1;
% figure(figcount)
% [pcx, pcy, pcz] = wakefilamentplot(1, ['C:\Users\hauti\heli-fw\MyCasesNew-' num2str(8) 'turns\TR\t' num2str(10) 'initial-mrtr\FWGEOM.DAT'],8);
% plot3(pcx(36,1:289,1)*rad(2),pcy(36,1:289,1)*rad(2),pcz(36,1:289,1)*rad(2),'LineWidth',1.5)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('10 kts - Case Ran with the MRTR Code')