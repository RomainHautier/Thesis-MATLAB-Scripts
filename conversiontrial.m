clear;
clc;
speed = 95;
file_combination = 'F';
turns = 8;

for i = 1:length(speed)

    trial = convert1to2rotors(speed(i),file_combination,turns);

if file_combination == 'H'
    fid = fopen(['C:\Users\hauti\heli-fwmrtr2\MyCases2rot\d' num2str(speed(i)) '\converted.txt'], 'wt');
    fprintf(fid, '%s', ['The wake files for ' num2str(speed(i)) 'kts have been converted on the ' datestr(now) '.' ...
        'This was realised using the 4 turn solution from HeliUM for the MR and 4 turn FW solution for the TR.']);
    fclose(fid);
elseif file_combination == 'F'
    fid = fopen(['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed(i)) '\converted.txt'], 'wt');
    fprintf(fid, '%s', ['The wake files for ' num2str(speed(i)) 'kts have been converted on the ' datestr(now) '.' ...
        ' This was realised by combining the solution for standalone runs of both the MR and TR with ' num2str(turns) ' Freewake turns.']);
    fclose(fid);
end

end