clear all, close all
d2r=pi/180;
speed = 10;
psi=0:10:360*d2r;
figcount = 1;

disp(['The circulation distribution over the rotor disks for at ' num2str(speed) ' kts.'])

% Loading the Standalone TR and MR data
FWG_b_TR = load(['C:\Users\hauti\heli-fwmrtr2\MyCasesNew\TR\t' num2str(speed) '\FWG_b.dat' ]);
FWG_b_MR = load(['C:\Users\hauti\heli-fwmrtr2\MyCasesNew\MR\t' num2str(speed) '\FWG_b.dat' ]);
%Loading the Standalone MR data from HeliUM
FWG_b_MRH = load(['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(speed) 'kt\FWG_B.dat']);

% % Loading and concantenating along the 3rd dimensions the 2 rotor cases - nn & nb 
% FWG_b_2Rnn = load(['C:\Users\hauti\heli-fwmrtr2\MyCases2rot\d' num2str(speed) '\nn-nzt577\FWG_b.dat' ]);
% FWG_b_2Rnb = load(['C:\Users\hauti\heli-fwmrtr2\MyCases2rot\d' num2str(speed) '\nb-nzt577\FWG_b.dat' ]);
% FWG_bs = cat(3,FWG_b_2Rnn, FWG_b_2Rnb);

%% Plotting the standalone TR circulation distribution

k=0;
for p = 1:36
    for s = 1:32
       k = k+1;
       LOCR(s,p) = FWG_b_TR(k,1);
       CIRCAFTER(s,p) = FWG_b_TR(k,2);
    end   
end

r = LOCR(:,1);
[th,r]=meshgrid(psi,r);
[x,y]=pol2cart(th,r);
CIRCAFTER(:,37) = CIRCAFTER(:,1);
figure(figcount)
contourf(x,y,CIRCAFTER)
colorbar;
title(['Circulation distribution over the standalone TR at ' num2str(speed) ' kts'])

figcount = figcount + 1;

%% Plotting the Standalone MR Circulation Distribution

k=0;
for p = 1:36
    for s = 1:32
       k = k+1;
       LOCR(s,p) = FWG_b_MR(k,1);
       CIRCAFTER(s,p) = FWG_b_MR(k,2);
    end   
end

r = LOCR(:,1);
[th,r]=meshgrid(psi,r);
[x,y]=pol2cart(th,r);
CIRCAFTER(:,37) = CIRCAFTER(:,1);
figure(figcount)
contourf(x,y,CIRCAFTER)
colorbar;
title(['Circulation distribution over the standalone MR at ' num2str(speed) ' kts'])

figcount = figcount + 1;

%% Plotting the circulation distribution over the MR & TR Rotors for the dissimilar rotor cases.

rotors = {'Main', 'Tail'};
conditions = {'nn', 'nb'};
for i = 1:length(FWG_bs(1,1,:))
    k=0;
    FWG_b = FWG_bs(:,:,i);

    % Arranging the data into the control points & circulation values.
    for r = 1:2
        for p = 1:36
            for s = 1:32
               k = k+1;
               LOCR(s,p) = FWG_b(k,1);
               CIRCAFTER(s,p) = FWG_b(k,2);
            end   
        end
        rad = LOCR(:,1);
        [th,rad]=meshgrid(psi,rad);
        [x,y]=pol2cart(th,rad);
        CIRCAFTER(:,37)=CIRCAFTER(:,1);
        figure(figcount)
        contourf(x,y,CIRCAFTER)
        colorbar
        figcount = figcount + 1;
        title([rotors{r} ' Rotor Circulation Distribution ' conditions{i}])
    end
end
