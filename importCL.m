function [CLr1,CLr2] = importCL(nr, speed)
% This function was built to import the Cl distribution in a single or 2
% rotor situation.
turns = 8;

if nr == 2
    filepath = ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed) '\nb\clhist.DAT']
else
    filepath = ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed) '\clhist.DAT']
end

mat = load([filepath]);
% mat(:,2:3) = [];

% Ordering the CL along each blade at each azimuthal position. The data for
% blade 1 is along the 1st layer of the z-dimension, blade 2 on the second
% layer and so on. A separate matrix is created for each rotor.

rows = 0;
seg_vec = 1:4:length(mat(:,1));

CLr1 = zeros(36,32,4);
CLr2 = CLr1;

rot1data_segmenting = [];
rot2data_segmenting = [];

for i = 2:2:length(seg_vec)-1

    rot1data_segmenting = [rot1data_segmenting; mat(seg_vec(i-1):seg_vec(i)-1,:)];
    rot2data_segmenting = [rot2data_segmenting; mat(seg_vec(i):seg_vec(i+1)-1,:)];

end

rot1data_segmenting(2877:2880,:) = mat(5753:5756,:);
rot2data_segmenting(2877:2880,:) = mat(5757:5760,:);


for j = 1:4
    index1 = find(rot1data_segmenting(:,2) == j);
    index2 = find(rot2data_segmenting(:,2) == j);
    CLr1(:,:,j) = rot1data_segmenting(index1(end-35:end,:),4:end);
    CLr2(:,:,j) = rot2data_segmenting(index1(end-35:end,:),4:end);
end

end