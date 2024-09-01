function [pcx, pcy, pcz] = wakefilamentplot(nr, filepath,turns)

np = 36; % Number of azimuthal steps
FW = load(filepath); % Loading the data
k = 0;

%% Check section: code to check if the concatenate of IWGEOM from the main and tail rotors
% into a single .data file is correct via visual representation.
% If check = 11, the check is performed on the wake of the two MR and TR
% standalone rotors.
% If check = 22, the check is performed on the wake of the two MR and TR
% wakes ran in the dissimilar rotor cases.

if nr == 2
    if turns == 4
        nzt = 145;
    elseif turns == 8
        nzt = 289;
    elseif turns == 16
        nzt = 577;
    end
    length(FW(:,1))/2;
    FW1 = FW(1:length(FW(:,1))/2,:);
    % disp(FW1(end,:))
    FW2 = FW(length(FW(:,1))/2+1:end,:);
    % disp(FW2(end,:))

    r = 1;
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
    
    k=0;
    r = 2;
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

%% Normal Plotting procedure
elseif nr == 1
    if turns == 4
        nzt = 145;
    elseif turns == 8
        nzt = 289;
    elseif turns == 16
        nzt = 577;
    end
    
    r = nr;
    for p=1:np
        for z=1:nzt
            k=k+1;
            psi(p,z,r)=FW(k,1);
            zeta(p,z,r)=FW(k,2);
            pcx(p,z,r)=FW(k,3);
            pcy(p,z,r)=FW(k,4);
            pcz(p,z,r)=FW(k,5);
        end
    end


end