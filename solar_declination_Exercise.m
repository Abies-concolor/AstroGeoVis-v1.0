%Solar declination Laboratory Exercise
%Instructor guidelines and suggestions are provided as code comments (lines beginning with percent sign)

%%Author: Dr. Tihomir Kostadinov, 2012 - 2020
%See accompanying Excel file for a simplified and modified version of this exercise

%Students are instructed to create an ascii file with their measurement as follows and load it here

%Columns of the input file loaded below should be arranged as follows: 
%Each column is a variable, and each row - an observation instance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Column 1): Year; 
%Column 2): Month number (1- Jan, 2 - Feb, etc.)
%Column 3): Day, e.g. enter 24 if the date is October 24th,
%Column 4): Observer latitude, decimal degrees, +ve in Northern Hemisphere,-ve in Southern Hemisphere
%Column 5): Observer longitude, decimal degrees, +ve in Eastern Hemisphere,-ve in Western Hemisphere
%Column 6): Observed shadow length, centimeters; 
%Column 7): Observed azimuth of the Sun, degrees (clockwise from true North); (Make sure to enter the azimuth of the direction of the Sun, NOT the direction the shadow is facing.)
%Column 8): True_Declination, degrees (looked up almanac value, see notes below) 
%Column 9): Length of the PVC pipe or other object used to create shadow, cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data courtesy of Fall 2012 GEOG 250 students, University of Richmond, Richmond, VA 
X = load('Solar_Declination_Exercise_Student_Data.txt'); %may need to modify path to file

%Notes: 
%1) Position of observer - latitude and longitude
%   Students can be given these coordinates or asked to determine them
%   themsleves, e.g. on Google Earth, or with a map or GPS unit
%   Note that sufficient accuracy is required.  A possible student task is to
%   investigate the propagation of error in these coordinates to the final
%   result. How much accuracy is sufficient? Note that longitude does not get used in the subsequent calculations; however, 
%   it is advisable for students to collect it for completeness. Students can be asked to discuss why longitude does not particpiate 
%   in the calculation of solar declination from solar altitude and azimuth.
%2) True solar declination - Students look up true calculated (predicted) declination, 
%   (e.g. on USNO Astronomical Applications or NOAA Solar Calculator websites),
%   or it can be calculated via the solar_coord.m script provided here) - this is the expected value that students compare with their
%   calculated value via the shadow length measurements. Best to use times
%   near the time of observation - students can be asked to discuss the time used to compute this valule as a source of uncertainty as well. . 
%3) To measure azimuth, students wil need to determine the direction of
%   true north sepaartely, or be given that information. 

%Instantaneous apparent altitude of the Sun at the moment of observation
%Students are either asked to figure out this formula in class/lab or as homework, or they are given the formula
shadow_length = X(:,6);
pipe_length = X(:,9);
h = atand(pipe_length./shadow_length); %degrees

%Computation of solar declination using conversin from horizontal to equatorial coordinates. 
%Students can be given this formula or asked to look it up from a reliable source.  Advanced math or astronomy students 
%can be asked to derive or prove it. Note that here azimuth is measured clockwise from true north, unlike in Meeus (1998)  

%See for example http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm (Positional Astronomy)
lat = X(:,4);
A = X(:,7);
sindelta = sind(h).*sind(lat) + (cosd(h).*cosd(lat)).*cosd(A); 
delta = asind(sindelta);

%Plot the observed vs almanac values of delta
delta_noaa = X(:,8);
figure
plot(delta_noaa,delta,'ro');
axis([-23.5 23.5 -23.5 23.5]);
hold on
plot([-23.5 23.5],[-23.5 23.5],'k--')
xticks([-20:5:20])
yticks([-20:5:20])
xlabel('Almanac Solar \delta, deg')
ylabel('Solar \delta from Shadow Length & Azimuth, deg')
title('Solar Declination Project Example Student Data')

%Students can also be asked to do/discuss the following:
%1) Add regression line to the plot and discuss performance metrics of their measurements, and outliers
% Should the regression be type I or type II (see #2 & #3 below)?
%2) What can be considered geophysical truth in an experiment? The almanac values come from models/compuations, NOT
% directly from measurements, so in some sense, their measurements can be considered as validation for almanac values.
% On the other hand, solar declination is an established value computed with very high precision by modern astronomical theory, 
% so their crude-instrument measurements are certain to contain higher
% error, by orders of magnitude.
%3) We often talk about declination on a given day, but declination is a
% continuously changing variable (i.e. it does not change in discrete steps).
% We are making an approximation if we consider a single value valid for
% the whole day, and error is intrioduced if the almanac value is not looked
% up for the exact time of the measurement. 

%4) Is declination following the expected trend for the time period covered by the project?

%5) Analyze possible sources of error, comment on regression stats, any outliers
% Some sources of error to consider (also see above): 
% 1) Measurements of the shadow length
    %pipe may not be vertical, ground may not be horizontal, actual
    %measurement itself (tape precision, the way you lay it out & other user errors)
% 2) Measureemnt of the azimuth - the protractor may not be oriented
    % perfectly E-W, actual angle measurement (protractor precisoion ,user errors) 

%Monte Carlo simulation for analysis and propagation of error: see accompanying script solar_declination_Monte_Carlo.m