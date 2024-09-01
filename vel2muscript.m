clear;
clc;

%% Rotor constants
global omTR RTR
omTR = 124.6; % TR angular velocity
RTR = 5.5; % Tail rotor radius
omMR = 27.00; % MR angular velocity
RMR = 26.83; % Main Rotor Radius (ft)
crdMR = 1.73; % Main rotor Chord (ft)
crdTR = 0.81; % Tail rotor chord (ft)
mMR = 86.7; % Main Rotor Blade Mass (slugs)

%% Rotation of velocities from the body frame to the TR frame + mu calcs.

% Setting up the rotation angles
a = deg2rad(0);
t = deg2rad(20);

% Body-frame TR velocity vector
us = [143.204	0	-2.9461

]';

% Setting up the transform matrix
% T = [1 0 0
%     0 cos(t) -sin(t)
%     0 sin(t) cos(t)];
T = [1 0 0
    0 1 0
    0 0 1];

t1 = deg2rad(70);
T1 = [1 0 0
    0 cos(90-t1) -sin(90-t1)
    0 sin(90-t1) cos(90-t1)];

% Calculating the velocities in the tail rotor reference frame.
ut = T*us
ut70 = T1*us

% Advance Ratio & Climb Advance Ratio at the TR.
mu = round(norm(ut(1:2),2)/(omTR*RTR),5);
disp(['mu = ' num2str(mu)])
muc = round(ut(3)/(omTR*RTR),5);
disp(['muc = ' num2str(muc)])

% Velocities in ft/sec
uvec = [168.76759	0	-2.125305]';

muMR = round(norm(uvec(1:2),2)/(omMR*RMR),5);
disp(['muMR = ' num2str(muMR)])
mucMR = round(uvec(3)/(omMR*RMR),5);
disp(['mucMR = ' num2str(mucMR)])

%% Main/Tail rotors inertia calculations
% Airfoil data is needed to calculate a density. We can get a mass per unit
% area

mpaMR = mMR/(crdTR*RMR) % Main rotor mass per unit area
mTR = mpaMR*crdTR*RTR

