clear;
clc;

BLKTR = 1; % Constant
BIASTR = 6; % Constant
S7 = 44645338.950671494; % Constant
S6 = 1.0188019125000000E-003;
rho = 2.0333962747827172E-003;
XMUXTR =  0.104391E+00;
XMUYTR = 0.319876E-02;
XMU2 = sqrt(XMUYTR^2+XMUXTR^2);
XMUZTR = 0.116425E-02;
V0TR = 0.0564357;
V0TR = 0.390362E-01;
XLAMTR = XMUZTR - V0TR;
VTTR = sqrt(XLAMTR^2+XMU2^2);
CT = 2*V0TR*VTTR;

%% Setting up the constants
target_CT = 0.0087;
TTR = CT*rho*S7*BLKTR;
theta = 0.3375; % Tail rotor blade pitch 
t0cTR = rad2deg(theta) + TTR*S6 - BIASTR

% while abs(CT-target_CT) >= 1e-4
% 
%     VTTR = target_CT/(2*V0TR);
% 
% 
% 
%     theta = 0.3347;
%     TTR = CT*rho*S7*BLKTR;
%     t0cTR = rad2deg(theta) + TTR*S6 - BIASTR
%     CT = 2*V0TR*VTTR;
% end
