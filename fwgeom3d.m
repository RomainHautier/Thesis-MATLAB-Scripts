clear;
close all;
% Plot geometry in 3D

set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontSize', 14);

%% Deciding whether to plot full helicopter dynamics rotor inflows or FreeWake ones.
% plotopt = 1: plot free wake rotor inflow for the main rotor,
% 2, plot the free wake rotor inflow for the tail rotor
%Define the velocity at which you want data to be taken.

% fwd_vel_vec = ['0'; '20'; '40'; '60'; '80'; '100']; 
fwd_vel_vec = 5:5:100;
fwd_vel_vec = [1 fwd_vel_vec];

% for figcount = 1:length(fwd_vel_vec)
%     fwd_vel = int2str(fwd_vel_vec(figcount));

% Load data and put into arrays

plotopt = 3; % Define which dataset to plot.
if plotopt == 1
    fn_s1 = 'C:\Users\hauti\heli-fw\MyCases\t'; % 'File name string 1'
elseif plotopt == 2
    fn_s1 = 'C:\Users\hauti\heli-fw\tmw-srt';
else
    fn_s1 = 'C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_';
end

fn_s2 = 'kt\FWGEOM.DAT'; % End of filename string.
% filename = convertCharsToStrings([fn_s1 fwd_vel fn_s2]);
filename = 'C:\Users\hauti\heli-fw2\MyCases2rot\d80\FWGEOM.DAT';
FW = load(filename);

np=36;
nzt=289;

% Here the position of each point along the 4 filaments studied are
% formatted into arrays.
% Each layer 'r' of the data represents the data obtained for one rotor.
% Each row 'p' of the data represents a timestep.
% Each column 'z' of the data represents one point of a filament.

%% Plotting the TWO IDENTICAL ROTOR CASES for both MFW and MRTR progs.
figcount = 80;
speed = 40:20:80;

for i = 1:length(speed)
    for j = 1:2
        figure(figcount)
        if j == 1
            fp1 = 'C:\Users\hauti\heli-fw2\MyCases2rot\d';
            fp2 = '\FWGEOM.DAT';
        else
            fp1 = 'C:\Users\hauti\heli-fwmrtr\MyCases2rotidentical\d';
            fp2 = 'notrim\FWGEOM.DAT';
        end

        [pcx, pcy, pcz] = wakefilamentplot(2,[fp1 num2str(speed(i)) fp2]);
        plot3(pcx(36,1:nzt,1),pcy(36,1:nzt,1),pcz(36,1:nzt,1),'LineWidth',1.5)
        hold on;
        plot3(pcx(36,1:nzt,2),pcy(36,1:nzt,2),pcz(36,1:nzt,2),'LineWidth',1.5)
        xlabel('x/R')
        ylabel('y/R')
        zlabel('z/R')
        figcount = figcount+1;

        if j == 1
            title(['FWD FLIGHT - MR & TR Wake Filaments at V = ' num2str(speed(i)) ' kts - MFW'])
        else
            title(['FWD FLIGHT - MR & TR Wake Filaments at V = ' num2str(speed(i)) ' kts - MRTR'])
        end
    end
end


%% GE

figcount = 100;
speed = [80 100];

for i = 1:length(speed)
    for j = 1:2
        figure(figcount)
        if j == 1
            fp1 = 'C:\Users\hauti\heli-fw2\MyCasesGE\d';
            fp2 = '\FWGEOM.DAT';
        else
            % fp1 = 'C:\Users\hauti\heli-fwmrtr\MyCases2GE\d';
            fp1 = 'C:\Users\hauti\heli-fwmrtr2\MyCasesGE\d';
            fp2 = '\FWGEOM.DAT';
        end

        [pcx, pcy, pcz] = wakefilamentplot(2,[fp1 num2str(speed(i)) fp2]);
        plot3(pcx(36,1:nzt,1),pcy(36,1:nzt,1),pcz(36,1:nzt,1),'LineWidth',1.5)
        hold on;
        plot3(pcx(36,1:nzt,2),pcy(36,1:nzt,2),pcz(36,1:nzt,2),'LineWidth',1.5)
        xlabel('x/R')
        ylabel('y/R')
        zlabel('z/R')
        figcount = figcount+1;

        if j == 1
            title(['GE - MR & TR Wake Filaments at V = ' num2str(speed(i)) ' kts - MFW'])
        else
            title(['GE - MR & TR Wake Filaments at V = ' num2str(speed(i)) ' kts - MRTR'])
        end
    end
end






%% MFW vs MRTR (trim vs no trim with mfw initial conditions)

filename1 = 'C:\Users\hauti\heli-fwmrtr\MyCases2rotidentical\d40notrim\FWGEOM.DAT';
FW1 = load(filename1);
filename2 = 'C:\Users\hauti\heli-fwmrtr\MyCases2rotidentical\d40notrim100ite\FWGEOM.DAT';
FW2 = load(filename2);

k=0;
for r = 1:2
    for p=1:np
        for z=1:nzt
            k=k+1;
            psi(p,z,r)=FW1(k,1);
            zeta(p,z,r)=FW1(k,2);
            pcx(p,z,r)=FW1(k,3);
            pcy(p,z,r)=FW1(k,4);
            pcz(p,z,r)=FW1(k,5);
        end
    end
end

figure(40)
plot3(pcx(1,1:nzt,1),pcy(1,1:nzt,1),pcz(1,1:nzt,1),'LineWidth',1.5)
hold on;
plot3(pcx(1,1:nzt,2),pcy(1,1:nzt,2),pcz(1,1:nzt,2),'LineWidth',1.5)
xlabel('x/R')
ylabel('y/R')
zlabel('z/R')
title(['Main & Tail Rotor Wake Filaments at V = ' filename1(38:40) ' kts'])
legend('Main Rotor', 'Tail Rotor')

k=0;
for r = 1:2
    for p=1:np
        for z=1:nzt
            k=k+1;
            psi(p,z,r)=FW2(k,1);
            zeta(p,z,r)=FW2(k,2);
            pcx(p,z,r)=FW2(k,3);
            pcy(p,z,r)=FW2(k,4);
            pcz(p,z,r)=FW2(k,5);
        end
    end
end

figure(41)
plot3(pcx(1,1:nzt,1),pcy(1,1:nzt,1),pcz(1,1:nzt,1),'LineWidth',1.5)
hold on;
plot3(pcx(1,1:nzt,2),pcy(1,1:nzt,2),pcz(1,1:nzt,2),'LineWidth',1.5)
xlabel('x/R')
ylabel('y/R')
zlabel('z/R')
title(['Main & Tail Rotor Wake Filaments at V = ' filename2(38:40) ' kts'])
legend('Main Rotor', 'Tail Rotor')


% for z=1:nzt
%     xx1(z)=pcx(1,z);
%     yy1(z)=pcy(1,z);
%     zz1(z)=pcz(1,z);
% end
% 
% % Fil 2
% 
% for z=1:nzt
%     xx2(z)=pcx(10,z);
%     yy2(z)=pcy(10,z);
%     zz2(z)=pcz(10,z);
% end
% 
% % Fil 3
% 
% for z=1:nzt
%     xx3(z)=pcx(19,z);
%     yy3(z)=pcy(19,z);
%     zz3(z)=pcz(19,z);
% end
% 
% % Fil 4
% 
% for z=1:nzt
%     xx4(z)=pcx(28,z);
%     yy4(z)=pcy(28,z);
%     zz4(z)=pcz(28,z);
% end
% 
% figure(figcount)
% plot3(xx1,yy1,zz1)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% plot3(xx2,yy2,zz2)
% hold on
% plot3(xx3,yy3,zz3)
% hold on
% plot3(xx4,yy4,zz4)
% 
% end
% 


