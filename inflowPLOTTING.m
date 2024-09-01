clear;
clc;
close all;

rotors = {'MR', 'TR'};
figcount = 1;

%% Calculating the inflow norm using the VBZ_PSI dat data.

% Defining the geometric constants required
omTR = 124.6;
omMR = 27;
RTR = 1.676;
RMR = 8.18;

% Defining logical constants
ns = 32;
np = 36;

% Array of speed over which the inflow needs to be extracted.
% speed = [10 15 20 40 50 60 80 95 100];
speed = [10:5:85 90 95 100];
% speed = 10:5:100;
% speed = [25 30 35 45 55 65 75 85];
saveDir = ['C:\Users\hauti\OneDrive - Imperial College London\Documents\London 23-24\Research Project\Result Figures\Inflow Plots\'];
%% Extracting the inflow from the VBZ_PSI.DAT file for both rotors MR & TR
% from the standalone FW cases. The first row of inflow belongs to the MR
% and the second row to the TR.
for i = 1:length(speed)

    %% Output matrix of the VBZPSI
    [VBZPSImatout, VBZPSImat, inflownorm] = importVBZPSI(1,speed(i),'');
    inflows_SA_VBZPSI(:,i) = inflownorm; % SA stands for 'standalone'

end 

%% Extracting the inflow of both rotors from VBZ_PSI.DAT from the 2 rotor cases

rotcase = {'nb'};
for case_index = 1:length(rotcase)
    for i = 1:length(speed)

        %% Output matrix of the VBZPSI
        [~, VBZPSImat, inflownorm] = importVBZPSI(2,speed(i),[rotcase{case_index} '\']);
        inflows_2rot_VBZPSI(:,i,case_index) = inflownorm;
    end 
end

% Fitting the inflow through a linear function for the 2 rotor case.
% Fitting a curve from 10 to 60 knots, one from 60 to 80 knots and 80 to
% 100 knots.

coeffs2rot = zeros(3,2);

fit_speed_limits = [10, 60, 85, 100];
for fs = 1:length(fit_speed_limits)
    inflow_ends(fs) = find(speed == fit_speed_limits(fs));
    speedfit{fs} = {0};
end

% Building a matrix with with the coefficients of the 3 fitting functions.
% Building vector of speeds over which to evaluate the fitting functions.
for in = 1:length(fit_speed_limits)-1
    % V = [];
    x = speed(inflow_ends(in):inflow_ends(in+1))';
    y = inflows_2rot_VBZPSI(2,inflow_ends(in):inflow_ends(in+1))';
    coeffs2rot(in,:) = polyfit(x,y,1);
    % 
    % % if in == 3
    % %     x0 = x(1);
    % %     y0 = y(1);
    % % else
    %     x0 = x(end);
    %     y0 = y(end);
    % % end
    % 
    % n = 1;
    % V(:,n+1) = ones(length(x),1,class(x));
    % for j = n:-1:1
    %     V(:,j) = x.*V(:,j+1);
    % end
    % A = [];
    % b = [];
    % Aeq = x0.^(n:-1:0);
    % coeffs2rot(in,:) = lsqlin(V,y,A,b,Aeq,y0);
    npoints = x(end) - x(1) + 1;
    % speedfit{in} = linspace(speed(inflow_ends(in))-2,speed(inflow_ends(in+1))+5,1000);
    speedfit{in} = speed(inflow_ends(in))-2:0.1:speed(inflow_ends(in+1))+5;
    fit{in} = polyval(coeffs2rot(in,:),speedfit{in});
end

%% Fitting the inflow through a linear function for the SA FW case.

coeffs2rotFW = zeros(3,2);

fit_speed_limits = [10, 40, 50, 90, 100];
for fs = 1:length(fit_speed_limits)
    inflow_ends(fs) = find(speed == fit_speed_limits(fs));
end

% Building a matrix with with the coefficients of the 3 fitting functions.
% Building vector of speeds over which to evaluate the fitting functions.
for in = 1:length(fit_speed_limits)-1
    V = [];
    x = speed(inflow_ends(in):inflow_ends(in+1))';
    y = inflows_SA_VBZPSI(2,inflow_ends(in):inflow_ends(in+1))';

    if in == 3
        x0 = x(1);
        y0 = y(1);
    else
        x0 = x(end);
        y0 = y(end);
    end

    n = 1;
    V(:,n+1) = ones(length(x),1,class(x));
    for j = n:-1:1
        V(:,j) = x.*V(:,j+1);
    end
    A = [];
    b = [];
    Aeq = x0.^(n:-1:0);
    coeffs2rotFW(in,:) = lsqlin(V,y,A,b,Aeq,y0);
    speedfitFW(in,:) = linspace(speed(inflow_ends(in))-2,speed(inflow_ends(in+1))+2,1000);
    fitFW(in,:) = polyval(coeffs2rotFW(in,:),speedfitFW(in,:));
end

% Doing the same for the FW standalone case, trying a 2nd degree fit. Here
% not end or starting point constraint is applied.
inflows_to_fit = {inflows_2rot_VBZPSI(2,:)',inflows_SA_VBZPSI(2,:)'};
for in = 1:2
    x = speed';
    y = inflows_to_fit{in};
    coeffs2nd(in,:) = polyfit(x,y,2);
    speedfit2(in,:) = linspace(10-2,100+2,1000);
    fit2(in,:) = polyval(coeffs2nd(in,:), speedfit2(in,:));
end


%% Extracting the MR & TR inflow from HeliUM (part output file and VBZ_PSI)

inflow_TR_helium = struct2array(load("inflow_TR_HeliUM.mat"));
speed_helium = [1 5:5:100];

for i = 1:length(speed_helium)
    %% Output matrix of the VBZPSI
    [~, VBZPSImat, inflownorm] = importVBZPSI(3,speed_helium(i),'');
    inflow_MR_helium(:,i) = inflownorm;
end 

%% Linear fitting of the TR inflow curve as per HeliUM.

coeffsTRHel = zeros(3,2);

fit_speed_limits = [10, 60, 80, 100];
for fs = 1:length(fit_speed_limits)
    inflow_ends(fs,:) = find(speed_helium == fit_speed_limits(fs));
end

% Building a matrix with with the coefficients of the 3 fitting functions.
% Building vector of speeds over which to evaluate the fitting functions.
for in = 1:length(fit_speed_limits)-1
    V = [];
    x = speed_helium(inflow_ends(in):inflow_ends(in+1))';
    y = inflow_TR_helium(inflow_ends(in):inflow_ends(in+1))';

    if in == 3
        x0 = x(1);
        y0 = y(1);
    else
        x0 = x(end);
        y0 = y(end);
    end

    n = 1;
    V(:,n+1) = ones(length(x),1,class(x));
    for j = n:-1:1
        V(:,j) = x.*V(:,j+1);
    end
    A = [];
    b = [];
    Aeq = x0.^(n:-1:0);
    coeffsTRHel(in,:) = lsqlin(V,y,A,b,Aeq,y0);
    speedfitHel(in,:) = linspace(speed_helium(inflow_ends(in))-2,speed_helium(inflow_ends(in+1))+2,1000);
    fitHel(in,:) = polyval(coeffsTRHel(in,:),speedfitHel(in,:));
end

x = speed_helium;
y = inflow_TR_helium;
coeffs2Hel = polyfit(x,y,2);
speedfit2hel = 10:0.1:100;
fitHel2 = polyval(coeffs2Hel, speedfit2hel);



%% Extracting the inflows of the FW cases from the timectcq.dat files.

% Standalone Main & Tail Rotors
for t = 1:length(rotors)
    [inflowmat, inflownorm] = importINFLOW(speed, ['C:\Users\hauti\heli-fwmrtr2\MyCasesNew-8turns\' rotors{t} '\t'],1,'');
    inflows_SA_timectcq(t,:) = inflownorm;
end


%% Plotting the MR inflows obtained from the VBZ_PSI.dat files.
figure(figcount)
plot(speed_helium,inflow_MR_helium,'color',[0 0.4470 0.7410],'LineWidth',1.5) % Helium Standalone MR
hold on;
scatter(speed_helium,inflow_MR_helium,18,[0 0.4470 0.7410],'filled') % Helium Standalone MR
hold on;
plot(speed,inflows_SA_VBZPSI(1,:),'color',[0.8500 0.3250 0.0980],'LineWidth',1.5) % FW Standalone MR
scatter(speed,inflows_SA_VBZPSI(1,:),18,[0.8500 0.3250 0.0980],'filled') % FW Standalone MR
for i = 1:length(inflows_2rot_VBZPSI(1,1,:))
    plot(speed,inflows_2rot_VBZPSI(1,:,i),'color',[0.4660 0.6740 0.1880],'LineWidth',1.5) % FW 2 Rot Cases
    scatter(speed,inflows_2rot_VBZPSI(1,:,i),18,[0.4660 0.6740 0.1880],'filled') % FW 2 Rot Cases
end
% plot(speed,inflows_SA_timectcq(1,:),'LineWidth',1.5) % FW Standalone MR - inflow(r,p) --> timectcq.dat

legend_vec_MR = {'HeliUM', 'FW - Standalone'};
for i = 1:length(rotcase)
    legend_vec_MR = [legend_vec_MR, ['trm = ' rotcase{i}(1) ', initial = ' rotcase{i}(2)]];
end
legend_vec_MR = [legend_vec_MR];

set(gca, 'YTick', 0:0.01:0.2)
set(gca, 'XTick', 0:10:100)
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
ylim([0 0.1])
grid on;

% title('Different Cases Inflow Plots','Interpreter','latex','FontSize',13)
xlabel('Forward Speed $V$ (kts)','Interpreter','latex','FontSize',14)
ylabel('Main Rotor Inflow $\lambda_{MR}$ (-)','Interpreter','latex','FontSize',14)
legend('HeliUM','', 'Freewake - Standalone Case','', 'Freewake - 2 Rotor Case (trm = ''n'', initial = ''b'')'  ,'Interpreter','latex','FontSize', 10)

figcount = figcount + 1;
%% Plotting the TR inflows obtained from the VBZ_PSI.dat files.
figure(figcount)
plot(speed_helium,inflow_TR_helium,'LineWidth',1.5) % Helium Standalone TR
hold on;
plot(speed,inflows_SA_VBZPSI(2,:),'LineWidth',1.5) % FW Standalone TR
for i = 1:length(inflows_2rot_VBZPSI(1,1,:))
    plot(speed,inflows_2rot_VBZPSI(2,:,i),'g-','LineWidth',1.5) % FW 2 Rot Cases
end
% plot(speed,inflows_SA_timectcq(2,:),'LineWidth',1.5) % FW Standalone MR - inflow(r,p) --> timectcq.dat

% Plotting the inflow fits.
for i = 1:length(fit)
    plot(speedfit{i},fit{i},'r--','Linewidth',1.3)
end

for i = 1:length(fitFW(:,1))
    plot(speedfitFW(i,:),fitFW(i,:),'k--','Linewidth',1.3)
end

% col = ['b','y'];
% for i = 1:2
%     plot(speedfit2(i,:), fit2(i,:),[col(i) '--'],'Linewidth',1.3)
% end

% for i = 1:3
%     plot(speedfitHel(i,:), fitHel(i,:), 'm--','Linewidth',1.3)
% end

plot(speedfit2hel,fitHel2,'c--','Linewidth',1.3)

xlabel('Forward Velocity (kts)')
ylabel('Average Inflow Norm (-) ')
title('Tail Rotor Inflow as a function of Forward Speed')
legend([legend_vec_MR 'Linear Fit - 2rot 1st deg','a','','','Quad Fit - 2rot','Linear Fit FW','Quad Fit FW'] ,'Location','southwest')

figcount = figcount + 1;


%% Finding the error between different fits to develop a speed dependent correction factor.

start_speed = 10;
end_speed = 100;
rot2fit_intercepts = [start_speed, 61, 88, end_speed];
SAfit_intercepts = [start_speed, 41, 50, 90.3, end_speed];

% Cleaning up the fit intervals where their segments intersect to obtain
% smooth transition between the different fit segments.
for i = 1:length(rot2fit_intercepts)-1
    rot2_error_intervals{i} = [rot2fit_intercepts(i):0.1:rot2fit_intercepts(i+1)-0.1; polyval(coeffs2rot(i,:),rot2fit_intercepts(i):0.1:rot2fit_intercepts(i+1)-0.1)];
end

for i = 1:length(SAfit_intercepts)-1
    SA_error_intervals{i} = [SAfit_intercepts(i):0.1:SAfit_intercepts(i+1)-0.1; polyval(coeffs2rotFW(i,:),SAfit_intercepts(i):0.1:SAfit_intercepts(i+1)-0.1)];
end



% Finding the error between the two
j1 = 1;
j2 = 1;
full_speed = 10:0.1:100;
counter = 1;
i1 = 1;
i2 = i1;
cond_hit_count = 0;
err_seg = zeros(1,1000);
for i = 1:length(10:0.1:100)
    full_speed(i)

    if SA_error_intervals{j1}(1,end) > rot2_error_intervals{j2}(1,end)
        if i2 > length(rot2_error_intervals{j2}(1,:))
            cond_hit_count = cond_hit_count+1;
            seg_length = length(find( err_seg(1:cond_hit_count-1,:) ~= 0));
            err_seg(cond_hit_count,1:length(rot2_error_intervals{j2}(1,:))) = factor_err(i-length(rot2_error_intervals{j2}(1,:)):i-1);
            % err_seg(cond_hit_count,i-seg_length) = factor_err(i-seg_length:i-1);
            speed_break(cond_hit_count) = full_speed(i); % Defining the speed at which the error calculated switches from one line fit to the next one.
            j2 = j2 + 1;
            i2 = 1;

        end
        err(i) = rot2_error_intervals{j2}(2,i2) - SA_error_intervals{j1}(2,i1);
        factor_err(i) = rot2_error_intervals{j2}(2,i2)/SA_error_intervals{j1}(2,i1);
        i1 = i1 + 1;
        i2 = i2 + 1;

    elseif i == length(full_speed)

        err(i) = rot2_error_intervals{j2}(2,end) - SA_error_intervals{j1}(2,end);
        factor_err(i) = rot2_error_intervals{j2}(2,end)/SA_error_intervals{j1}(2,end);
        err_seg(cond_hit_count,1:length(SA_error_intervals{j1}(1,:))) = factor_err(i-length(SA_error_intervals{j1}(1,:)):i-1);

    else
        if i1 > length(SA_error_intervals{j1}(1,:))
            cond_hit_count = cond_hit_count+1;
            % err_seg(cond_hit_count,1:length(SA_error_intervals{j1}(1,:))) = err(i-length(SA_error_intervals{j1}(1,:)):i-1);
            err_seg(cond_hit_count,1:length(SA_error_intervals{j1}(1,:))) = factor_err(i-length(SA_error_intervals{j1}(1,:)):i-1);
            speed_break(cond_hit_count) = full_speed(i); % Defining the speed at which the error calculated switches from one line fit to the next one.
            j1 = j1 + 1;
            i1 = 1;
         
        end
        err(i) = rot2_error_intervals{j2}(2,i2) - SA_error_intervals{j1}(2,i1);
        factor_err(i) = rot2_error_intervals{j2}(2,i2)/SA_error_intervals{j1}(2,i1);
        i1 = i1 + 1;
        i2 = i2 + 1;
    end
    counter = counter + 1;
end

int1 = [];
int2 = [];
for i = 1:length(rot2_error_intervals)
    int1 = [int1 rot2_error_intervals{i}(2,:)];
end

for i = 1:length(SA_error_intervals)
    int2 = [int2 SA_error_intervals{i}(2,:)];
end

facts = int1./int2;% Correction factor tabulated from SA to 2 rotor case.
facts(901) = facts(end);
int2(901) = int2(end);
facts2 = int2./fitHel2;

figure(figcount)

tiles = tiledlayout("vertical");
xlabel(tiles,'Forward Speed $V$ (kts)','Interpreter','latex','FontSize',16)
ylabel(tiles,'Tail Rotor Inflow $\lambda_{TR}$ (-)','Interpreter','latex','FontSize',16)

% First Tile - HeliUM TR Inflowss

nexttile
plot(speed_helium,inflow_TR_helium,'color',[0 0.4470 0.7410],'LineWidth',1.5) % Helium Standalone TR
hold on;
% Plotting the HeliUM TR inflow fit.
plot(speedfit2hel,fitHel2,'color',[0 0.4470 0.7410],'Linewidth',1.6,'LineStyle','--')
legend('$\lambda_{0}$','Polyfit of $\lambda_{0}$','Interpreter','latex','FontSize', 10)
title('HeliUM','Interpreter','latex','FontSize',13)
ylim([0 0.1])

% Second tile - SA TR Inflows
nexttile
plot(speed,inflows_SA_VBZPSI(2,:),'color',[0.8500 0.3250 0.0980],'LineWidth',1.5) % FW Standalone TR
hold on;
for i = 1:length(SA_error_intervals)
    plot(SA_error_intervals{i}(1,:),SA_error_intervals{i}(2,:),'color',[0.8500 0.3250 0.0980],'LineWidth',1.6,'LineStyle','--')
end
legend('$\lambda_{avg}$','Polyfit of $\lambda_{avg}$','Interpreter','latex','FontSize', 10)
title('Freewake - Standalone Rotor Case','Interpreter','latex','FontSize',13)
ylim([0 0.1])

% Third Tile - 2 Rotor TR Inflows
nexttile

for i = 1:length(inflows_2rot_VBZPSI(1,1,:))
    plot(speed,inflows_2rot_VBZPSI(2,:,i),'color',[0.4660 0.6740 0.1880],'LineWidth',1.5) % FW 2 Rot Cases
    hold on;
end
for i = 1:length(rot2_error_intervals)
    plot(rot2_error_intervals{i}(1,:),rot2_error_intervals{i}(2,:),'color',[0.4660 0.6740 0.1880],'LineWidth',1.6,'LineStyle','--')
end
legend('$\lambda_{avg}$','Polyfit of $\lambda_{avg}$','Interpreter','latex','FontSize', 10)
title('Freewake - Two Rotor Case','Interpreter','latex','FontSize',13)
ylim([0 0.2])

allAxes = findall(gcf,'Type','axes'); % Find all axes in the figure
set(allAxes, 'XLim', [10 100]);         % Set common xlim
set(allAxes, 'XGrid', 'on', 'YGrid', 'on');
set(allAxes, 'XTick', 10:10:100) %sets the x-tick marks at intervals of 2 from 0 to 10 for all tiles.
set(allAxes, 'XMinorTick', 'on')
set(allAxes, 'YMinorTick', 'on') %enable minor ticks for both x and y axes.

figcount = figcount + 1;

% Proper Plotting of all inflows.

figure(figcount)
plot(speed_helium,inflow_TR_helium,'color',[0 0.4470 0.7410],'LineWidth',1.5) % Helium Standalone TR
hold on;
scatter(speed_helium,inflow_TR_helium,18,[0 0.4470 0.7410],'filled')

hold on;
plot(speed,inflows_SA_VBZPSI(2,:),'color',[0.8500 0.3250 0.0980],'LineWidth',1.5) % FW Standalone TR
scatter(speed,inflows_SA_VBZPSI(2,:),18,[0.8500 0.3250 0.0980],'filled')
hel_vs_fw_error = mean(abs(inflow_TR_helium(3:end) - inflows_SA_VBZPSI(2,:)'));
hel_vs_fw_error_percentage = mean(abs(inflow_TR_helium(3:end) - inflows_SA_VBZPSI(2,:)')./inflow_TR_helium(3:end));

inflows2rot = [];
for i = 1:length(inflows_2rot_VBZPSI(1,1,:))
    plot(speed,inflows_2rot_VBZPSI(2,:,i),'color',[0.4660 0.6740 0.1880],'LineWidth',1.5) % FW 2 Rot Cases
    inflows2rot = [inflows2rot inflows_2rot_VBZPSI(2,:,i)];
    scatter(speed,inflows_2rot_VBZPSI(2,:,i),18,[0.4660 0.6740 0.1880],'filled') % FW 2 Rot Cases
    hold on;
end

hel_vs_2rot_error = mean(abs(inflow_TR_helium(3:end) - inflows2rot'));
hel_vs_2rot_error_percentage = mean(abs(inflow_TR_helium(3:end) - inflows2rot')./inflow_TR_helium(3:end));
set(gca, 'YTick', 0:0.04:0.2)
set(gca, 'XTick', 0:10:100)
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
grid on;
legend('HeliUM','', 'Freewake - Standalone Case','', 'Freewake - 2 Rotor Case (trm = ''n'', initial = ''b'')'  ,'Interpreter','latex','FontSize', 10)
% title('Different Cases Inflow Plots','Interpreter','latex','FontSize',13)
xlabel('Forward Speed $V$ (kts)','Interpreter','latex','FontSize',14)
ylabel('Tail Rotor Inflow Norm $\lambda_{avg}$ (-)','Interpreter','latex','FontSize',14)
ylim([0 0.16])

figcount = figcount + 1;


%% Writing the file to be read by HeliUM to introduce the inflow correction.

vel = round([10:0.1:100]',2);
correct_SA_2rot = int2.*facts;
correct_Hel_SA = fitHel2.*facts2;
full_correction = fitHel2.*facts2.*facts;
tabulated_data = [vel, facts2', facts'];

figure(figcount)
plot(vel,full_correction,'LineWidth',1.6)
hold on;
for i = 1:length(inflows_2rot_VBZPSI(1,1,:))
    plot(speed,inflows_2rot_VBZPSI(2,:,i),'color',[0.4660 0.6740 0.1880],'LineWidth',1.5) % FW 2 Rot Cases
    scatter(speed,inflows_2rot_VBZPSI(2,:,i),18,[0.4660 0.6740 0.1880],'filled') % FW 2 Rot Cases
    hold on;
end
title('Verification of the full correction')




% Open a file for writing
fileID = fopen('inflowcorrectiondata.DATA', 'w');

for i = 1:length(vel)
    fprintf(fileID, '%.1f %f %f\n', tabulated_data(i,1), tabulated_data(i,2), tabulated_data(i,3));
end

% fileID = fopen('inflowcorrectiondata.DATA','w');
% for i = 1:length(vel)
%     fprintf(fileID, '%0.1f %f\n', [vel(i) full_correction(i)]);
% end


% % Write the data to the file
% % %f is for floating-point numbers, separated by spaces
% fprintf(fileID, '%f %f %f\n', tabulated_data);

% Close the file
fclose(fileID);

% plot(10:0.1:100, correct_Hel_SA,"Color",[0.2 0.1 0.4],'LineWidth',1.3,'LineStyle','--')
% plot(10:0.1:100, correct_SA_2rot,"Color",[0.2 0.6 0.4],'LineWidth',1.3,'LineStyle','--')
% plot(full_speed,err,'r--','LineWidth',1.6)
% legend('HeliUM TR','FW SA','2 Rotor Case -','HeliUM TR Fit','SA Tail Rotor Fit','','','2 Rotor Inflow Fit' ,'','','','HeliUM to SA Corrected Fit', 'SA to 2 Rotor Corrected fit','Error Between Fits')

% figcount = figcount + 1;


% figure(8)
% plot(10:0.1:99.9,facts)


% %% Linear Fitting the different segments of error fitting
% speed_break = [10 speed_break 100.1];
% sum_index = 1;
% err_seg = {};
% for s = 1:length(speed_break)-1
%     x = speed_break(s):0.1:speed_break(s+1)-0.1;
%     err_seg{s} = factor_err(sum_index:sum_index+length(x)-1);
%     FW2rot_factor_fun(s,:) = polyfit(x,err_seg{s},1);
%     sum_index = sum_index + length(x)-1;
%     xfit{s} = x;
%     FW2rot_corrected{s} = polyval(FW2rot_factor_fun(s,:),xfit{s});
%     plot(xfit{s},FW2rot_corrected{s}.*SA_error_intervals{1}(2,:),'LineWidth',1.6,'Color','y')
% end
% 
% figure(figcount)
% plot(speed_helium,inflow_TR_helium,'LineWidth',1.5) % Helium Standalone MR
% hold on;
% plot(speed,inflows_SA_VBZPSI(2,:),'LineWidth',1.5) % FW Standalone TR
% for i = 1:length(inflows_2rot_VBZPSI(1,1,:))
%     plot(speed,inflows_2rot_VBZPSI(2,:,i),'g-','LineWidth',1.5) % FW 2 Rot Cases
% end
% 
% % Plotting the inflow fits.
% for i = 1:length(fit)
%     plot(speedfit{i},fit{i},'r--','Linewidth',1.3)
% end
% 
% for i = 1:length(fitFW(:,1))
%     plot(speedfitFW(i,:),fitFW(i,:),'k--','Linewidth',1.3)
% end
% 
% for i = 1:length(FW2rot_corrected)
%     plot(xfit{i},FW2rot_corrected{i},'c--','Linewidth',1.3)
% end






%% Applying the same methodology to apply a correction factor from Hel SA to FW SA








