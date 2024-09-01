function [us, mu, muc, ut, muMR, mucMR] = vel2mu(TR_velb,hel_vel)
%% Rotor constants
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


% Setting up the transform matrix
T = [1 0 0
    0 cos(t) -sin(t)
    0 sin(t) cos(t)];

for i = 1:length(TR_velb(1,:))
    
    us(:,i) = TR_velb(:,i);
    % Calculating the velocities in the tail rotor reference frame.
    ut(i,:) = T*us(:,i);
    
    % Advance Ratio & Climb Advance Ratio at the TR.
    mu(i) = round(norm(ut(i,1:2),2)/(omTR*RTR),5);
    muc(i) = round(ut(i,3)/(omTR*RTR),5);

    % Advance Ratio & Climb Advance Ratio at the MR.
    uvec = hel_vel(:,i);
    muMR(i) = round(norm(uvec(1:2),2)/(omMR*RMR),5);
    mucMR(i) = round(uvec(3)/(omMR*RMR),5);


end



end