clear;
% clc;
close all;

filepath = 'C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\8turns';
input_speed = 10;
output_speed = 20;
solution = 'nb';
state = DATtoDATAconverter(filepath, num2str(input_speed), num2str(output_speed), solution);