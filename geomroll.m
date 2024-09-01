% Tilt rear geometry data for a roll maneuver

function [pcynew,pcznew]=geomroll(pcy,pcz,roll)

n=length(pcy);

roll=roll*pi/180;

for i=1:n
    pcynew(i)=pcy(i)*cos(roll)+pcz(i)*sin(roll);
    pcznew(i)=-pcy(i)*sin(roll)+pcz(i)*cos(roll);
end

% NOTE:
%
% The roll angle is entered in body axes, that is, positive to the right.
% The geometry is in wake coordinate system, with x pointing back and
% therefore roll is positive to the left. This is already taken into
% account in the transformation.