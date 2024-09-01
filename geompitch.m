% Tilt rear geometry data for a roll maneuver

function [pcxnew,pcznew]=geompitch(pcx,pcz,pitch)

n=length(pcx);

pitch=pitch*pi/180;

for i=1:n
    pcxnew(i)=pcx(i)*cos(pitch)+pcz(i)*sin(pitch);
    pcznew(i)=-pcx(i)*sin(pitch)+pcz(i)*cos(pitch);
end

