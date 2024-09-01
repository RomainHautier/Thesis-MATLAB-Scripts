% Plot the geometry, with the blades on, in the three views
% with the actual body attitudes

function plotgeomangles(roll,pitch)

close all %, clear all
set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontSize', 14);

% Orienting the wake wrt roll and pitch angles.
if nargin==0
    roll=0;
    pitch=0;
end

%% Deciding whether to plot full helicopter dynamics rotor inflows or FreeWake ones.
% plotopt = 1: plot free wake rotor inflow for the main rotor,
% 2, plot the free wake rotor inflow for the tail rotor

%Define the velocity at which you want data to be taken.
fwd_vel = '100'; 

% Load data and put into arrays

plotopt = 3; % Define which dataset to plot.
if plotopt == 1
    fn_s1 = 'C:\Users\hauti\heli-fw\tmw-srm'; % 'File name string 1'
elseif plotopt == 2
    fn_s1 = 'C:\Users\hauti\heli-fw\tmw-srt';
else
    fn_s1 = 'C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_';
end

fn_s2 = ['kt\FWS.DAT';'kt\FWR.DAT';'kt\FWT.DAT']; % End of filename string.
% Plotting options depending on the fwd velocity of the helicopter.

for i = 1:3
    filename(i,:) = convertCharsToStrings([fn_s1 fwd_vel fn_s2(i,:)]);
end

side=load(filename(1,:));
rear=load(filename(2,:));
top = load(filename(3,:));

% side=load('C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_100kt\FWS.DAT');
% rear=load('C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_100kt\FWR.DAT');
% top = load('C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_100kt\FWT.DAT');

pcx=side(:,1);
pcy=rear(:,1);
pcz=rear(:,2);
 
[pcy1,pcz1]=geomroll(pcy,pcz,roll);      % REAR
[pcx2,pcz2]=geompitch(pcx,pcz,pitch);    % SIDE

% NOTE:
%
% The roll angle is entered in body axes, that is, positive to the right.
% The geometry is in wake coordinate system, with x pointing back and
% therefore roll is positive to the left. This is already taken into
% account in the transformation.
   
   % Define blades

bladestop1=[0,0;pcx(1),pcy(1)];
bladestop2=[0,0;pcx(146),pcy(146)];
bladestop3=[0,0;pcx(291),pcy(291)];
bladestop4=[0,0;pcx(436),pcy(436)];

bladesside1=[0,0;pcx2(1),pcz2(1)];
bladesside2=[0,0;pcx2(146),pcz2(146)];
bladesside3=[0,0;pcx2(291),pcz2(291)];
bladesside4=[0,0;pcx2(436),pcz2(436)];

bladesrear1=[0,0;pcy1(1),pcz1(1)];
bladesrear2=[0,0;pcy1(146),pcz1(146)];
bladesrear3=[0,0;pcy1(291),pcz1(291)];
bladesrear4=[0,0;pcy1(436),pcz1(436)];


% TOP VIEW

figure
plot(pcx(1:145),pcy(1:145),'-k','LineWidth',1)
hold on
plot(pcx(146:290),pcy(146:290),'-k','LineWidth',1)
hold on
plot(pcx(291:435),pcy(291:435),'-k','LineWidth',1)
hold on
plot(pcx(436:580),pcy(436:580),'-k','LineWidth',1)
hold on
plot(bladestop1(:,1),bladestop1(:,2),'-k','LineWidth',3)
hold on
plot(bladestop2(:,1),bladestop2(:,2),'-k','LineWidth',3)
hold on
plot(bladestop3(:,1),bladestop3(:,2),'-k','LineWidth',3)
hold on
plot(bladestop4(:,1),bladestop4(:,2),'-k','LineWidth',3)
axis([-1 3 -1.5 1.5])
%axis([-1 3 -2 2])
xlabel('x/R','FontName', 'Helvetica','FontSize', 14)
ylabel('y/R','FontName', 'Helvetica','FontSize', 14)
%set(gca,'XTick',[-1,0,1,2,3])
%set(gca,'YTick',[-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
%set(gca,'YTick',[-2,-1,0,1,2])
%set(gca,'YTickLabel',['-1.5|-1.0|-0.5|0.0|0.5|1.0|1.5'])
%set(gca,'YTickLabel',['-1.5','-1.0','-0.5','0.0','0.5','1.0','1.5'])
print -depsc topview.eps

% SIDE VIEW

figure
plot(pcx2(1:145),pcz2(1:145),'-k','LineWidth',1)
hold on
plot(pcx2(146:290),pcz2(146:290),'-k','LineWidth',1)
hold on
plot(pcx2(291:435),pcz2(291:435),'-k','LineWidth',1)
hold on
plot(pcx2(436:580),pcz2(436:580),'-k','LineWidth',1)
hold on
plot(bladesside1(:,1),bladesside1(:,2),'-k','LineWidth',3)
hold on
plot(bladesside2(:,1),bladesside2(:,2),'-k','LineWidth',3)
hold on
plot(bladesside3(:,1),bladesside3(:,2),'-k','LineWidth',3)
hold on
plot(bladesside4(:,1),bladesside4(:,2),'-k','LineWidth',3)
%axis([-1 3 -1 1])
axis([-1 3 -1.5 1.5])
xlabel('x/R','FontName', 'Helvetica','FontSize', 14)
ylabel('z/R','FontName', 'Helvetica','FontSize', 14)
%set(gca,'XTick',[-1,0,1,2,3])
%set(gca,'YTick',[-1.0,-0.5,0,0.5,1.0])
%set(gca,'YTickLabel',['-1.0|-0.5|0.0|0.5|1.0'])
%set(gca,'YTick',[-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
%set(gca,'YTickLabel',['-1.5|-1.0|-0.5|0.0|0.5|1.0|1.5'])
print -depsc sideview.eps

% REAR VIEW

figure
plot(pcy1(1:145),pcz1(1:145),'-k','LineWidth',1)
hold on
plot(pcy1(146:290),pcz1(146:290),'-k','LineWidth',1)
hold on
plot(pcy1(291:435),pcz1(291:435),'-k','LineWidth',1)
hold on
plot(pcy1(436:580),pcz1(436:580),'-k','LineWidth',1)
hold on
plot(bladesrear1(:,1),bladesrear1(:,2),'-k','LineWidth',3)
hold on
plot(bladesrear2(:,1),bladesrear2(:,2),'-k','LineWidth',3)
hold on
plot(bladesrear3(:,1),bladesrear3(:,2),'-k','LineWidth',3)
hold on
plot(bladesrear4(:,1),bladesrear4(:,2),'-k','LineWidth',3)
%axis([-1.5 1.5 -1 1])
axis([-1.5 1.5 -1.5 1.5])
xlabel('y/R','FontName', 'Helvetica','FontSize', 14)
ylabel('z/R','FontName', 'Helvetica','FontSize', 14)
%set(gca,'XTick',[-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
%set(gca,'XTickLabel',['-1.5|-1.0|-0.5|0.0|0.5|1.0|1.5'])
%set(gca,'YTick',[-1.0,-0.5,0,0.5,1.0])
%set(gca,'YTickLabel',['-1.0|-0.5|0.0|0.5|1.0'])
%set(gca,'YTick',[-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
%set(gca,'YTickLabel',['-1.5|-1.0|-0.5|0.0|0.5|1.0|1.5'])
print -depsc rearview.eps