clear;
close all;
clc;

%% Setting up constants

omTR = 124.6; % TR angular velocity
RTR = 5.5; % Tail rotor radius
omMR = 27.00; % MR angular velocity
RMR = 26.83; % Main Rotor Radius (ft)
crdMR = 1.73; % Main rotor Chord (ft)
crdTR = 0.81; % Tail rotor chord (ft)
mMR = 86.7; % Main Rotor Blade Mass (slugs)
asr=20*pi/180;
m2f=3.28084; % Meters to feet conversion
ns = 32; % Number of blade segments
np = 36; % Number of timesteps

%% Setting the range of speeds over which to examine the results.
speed_helium = [1 5:5:100];
trimdata = importTRIMDATA(speed_helium);


%% Importing Excel data table from HeliUM runs output files.
data = importEXCEL("C:\Users\hauti\MATLAB Drive\Research Project\FWHELsheet.xlsx");
CTs = data(:,2:3);
% CT_MR = CTs(:,1);
CT_TR = flip(CTs(:,2));
% hel_vel = data(:,4:6)';
% P_MR = data(:,7);
% T_MR = data(:,8);
lambda_TR = data(:,9);
% uc = data(:,10:11);
% uc_MR = uc(:,1);
% uc_TR = uc(:,2);
t0r = data(:,12); % This output comes from the AATAILR.DEBUG file
TR_velb = data(:,13:15)';
t0rnotrim = data(:,21);
TTR = data(:,22);

%% Retrieving data from the TRIMPLOT Data file.
vel = trimdata(:,1);
c_angle = trimdata(:,2);
turn_rate = trimdata(:,3);
P_MR = trimdata(:,4);
T_MR = trimdata(:,5);
CT_MR = trimdata(:,6);
CT_sigma = trimdata(:,7);

% The next 4 variables are control inputs in inches, need to convert them
% to degrees using the UH-60a control-input conversion
t1c = trimdata(:,8);
t1s = trimdata(:,9);
t0r_MR = trimdata(:,10);
t0r_helium_ouput = trimdata(:,11); % Tail rotor collective from TRIMPLOTDATA.FLX
hel_vel = trimdata(:,12:14)';

%% Converting the control inputs from inches to degrees

% Main Rotor Collective
neutral_t0r = 9.9; % Neutral collective pitch angle of the MR (deg)
t0r_MR = neutral_t0r + t0r_MR.*16/10;

% Tail Rotor Collective
BIASTR = 6; % Bias angle on collective pitch (deg)
DELTTR = 0.145500e-02; % Rate of change of coning with thrust.
TD3TR = 0.700207; % Flapping hinge offset angle (deg)
S6 = DELTTR*TD3TR;
neutral_t0r_TR_deg = 15; % Neutral collective pitch angle of the TR (deg)
neutral_t0r_TR_in = 2.69; % Neutral collective pitch angle of the TR (inches)
t0r_TR_pedal = neutral_t0r_TR_deg + (t0r_helium_ouput - neutral_t0r_TR_in).*29.8/5.38; % Tail rotor collective pitch input conversion from inches to degrees.
t0r_TR = t0r_TR_pedal + flip(TTR.*S6 + BIASTR);


% Main rotor lateral cyclic
neutral_t1c_deg = -1.54;
neutral_t1c_in = 5.96;

% Main rotor longitudinal cyclic
neutral_t1s_deg = -11.0;
neutral_t1s_in = 8.89;

for i = 1:length(t1c)
    if t1c(i) <= neutral_t1c_in
        t1c(i) = -8.00 + t1c(i)*(8+neutral_t1c_deg)/5;
    else
        t1c(i) = neutral_t1c_deg + (t1c(i)-5)*(8-neutral_t1c_deg)/5;
    end

    if t1s <= neutral_t1s_in
        t1s(i) = 16.50 - t1s(i)*(16.5-neutral_t1s_deg)/5;
    else
        t1s(i) = neutral_t1s_deg - (t1s(i)-5)*(12.3+neutral_t1s_deg)/5;
    end
end

[us, mu, muc, TR_vel, muMR, mucMR] = vel2mu(TR_velb,hel_vel);
recap = table(vel, flip(mu'), flip(muc'), muMR', mucMR', t0r_MR, t1c, t1s, t0r_TR,'VariableNames', {'Velocity','muTR'  'mucTR'  'muMR'  'mucMR' 'MR Collective (deg)'  'MR Lateral Cyclic (deg)' 'MR Longitudinal Cyclic (deg)' 'TR Collective (deg)'});
disp(recap)


figure_count = 1;

figure(figure_count)
plot(vel, P_MR,'LineWidth',2)
xlabel('Forward Velocity (kts)')
ylabel('Main Rotor Power (HP)')

figure_count = figure_count+1;

figure(2)
plot(vel, CT_TR/max(CT_TR),'LineWidth',2)
hold on;
plot(vel, CT_MR/max(CT_MR),'LineWidth',2)
xlabel('Forward Velocity (kts)')
ylabel('Normalised Rotor Thrust Coefficient (-)')
legend('Tail Rotor','Main Rotor')
figure_count = figure_count+1;

figure(3)
plot(flip(vel),lambda_TR,'LineWidth',2)
xlabel('Forward Velocity (kts)')
ylabel('Tail Rotor Inflow Coefficient (-)')
figure_count = figure_count+1;

figure(4)
plot(vel, T_MR,'LineWidth',2)
xlabel('Forward Velocity (kts)')
ylabel('Main Rotor Thrust (lbs)')
figure_count = figure_count+1;

%% Extracting data from the Freewake files
speed = 10:5:100;
% Extracting the inflow information from the timectq.dat file and
% extracting the norm of the inflow at each velocity.
% [inflowmat, inflownorm] = importINFLOW(speed, 'C:\Users\hauti\heli-fw\MyCasesNew\t',1,'');
% [inflowmatnot1cs, inflownormnot1cs] = importINFLOW(speed, 'C:\Users\hauti\heli-fw\MyCasesnot1cs\t',1,'');
% Extracting information from the timelfightmodified.dat file to obtain the
% tail rotor collective pitch angle when it is ran on its own.

% [timeflightmat, timeflight, t0rfw, t1cfw, t1sfw, b0rTR, b1cTR, b1sTR] = importtimeflight(speed(2:end),'C:\Users\hauti\heli-fw\MyCasesnot1cs\t','');

speed = 10;
for i = 1:length(speed)

    %% Output matrix of the VBZPSI
    VBZPSImat = importVBZPSI(speed(i));
    VBZpsi(1:length(VBZPSImat(:,1)),1:length(VBZPSImat(1,:)),i) = VBZPSImat;

end 

for i = 1:length(speed)

    VBZ(i,:) = VBZpsi(:,6,i);
    
    s = 32:32:length(VBZ(1,:));
    k=1;
    
    for p = 1:np
        inflownw(1:ns,p,i) = VBZ(i,k:s(p));
        k = k+ns;
    end
    norm(i) = sum(sum(inflownw(:,:,i)));

end

normmps=-norm/(ns*np);
normfps=normmps*m2f;
lambda_TRfw=normmps/(omTR*RTR);

figure(figure_count)
plot(flip(vel),t0r,'LineWidth',2)
hold on;
% plot(vel(3:end),t0rfw,'LineWidth',2)
plot(flip(vel(1:length(t0rnotrim))),t0rnotrim,'LineWidth',2)
xlabel('Forward Velocity (kts)')
ylabel('Tail Rotor Collective Pitch Angle (deg)')
legend('HeliUM','Freewake (standalone rotor)','Freewake no trim (t1c = t1s = 0)')
figure_count = figure_count+1;

figure(figure_count)
plot(flip(vel), lambda_TR, 'LineWidth',2)
hold on;
plot(vel(3:end),lambda_TRfw,'LineWidth',2)
plot(vel(3:end),inflownorm,'LineWidth',2)
plot(vel(3:end),inflownormnot1cs,'LineWidth',2)
xlabel('Forward Velocity (kts)')
ylabel('Tail Rotor Inflow Coefficient (-)')
legend('HeliUM','Freewake (standalone rotor)- calculated from VBZPSI','Freewake (standalone rotor)- extracted from inflow(r,p) written in timectcq.dat','Freewake no trim (t1c = t1s = 0)')
figure_count = figure_count+1;

figure(figure_count)
plot(vel,t0r_MR,'LineWidth',2)
hold on;
plot(vel, t0r_TR)
