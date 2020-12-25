function [s, So] = sidereal_time(JD,long)
% Based on the algorithms of Meeus (1998), Astronomical Algorithms, Ch. 12, Ch. 13, Ch. 22
%
% INPUT: 
% JD: Julian date number, can be a vector
% long: geographic longitude in degrees, positive in Eastern Hemisphere
%       Note that this convention differs from the one adopted in Meeus (1998)  - see his Ch. 13
%       long should be a scalar
%
% OUTPUT:
% s: the apparent local sidereal time in HOURS
% So: the apparent sidereal time for the given instant at Greenwich, in hours
%
% Caution! Input not checked for validity in this function. 
%
% Implemented in MATLAB(r) and vectorized by Dr. Tihomir S. Kostadinov, 2014 - 2020
% Substantial changes implemented in Fall 2020

T = (JD - 2451545.0)./36525; 

%Mean sidereal time at Greenwich - Meeus Eq. 12.4
So = 280.46061837 + 360.98564736629*(JD-2451545.0) + 0.000387933*T.^2 - (T.^3)/38710000; %in degrees

%Correcting for nutation, see pg. 88, Ch. 12 of Meeus (1998) 
omega = 125.04452 - 1934.136261.*T; %degrees, Ch. 22, pg 144

%L = 280.4665 + 36000.7698*T; %Mean longitude of the Sun, degrees
L = 280.46646 +36000.76983*T + 0.0003032*(T.^2); %geometric mean longitude of the Sun, referred to the mean equinox of date, as given in Ch. 25
L_prime = 218.3165 + 481267.8813*T; %Mean longitude of the Moon, degrees

%nutation in longitude
delta_psi_arcsec = -17.20*sind(omega) - 1.32*sind(2*L) - 0.23*sind(2*L_prime) + 0.21*sind(2*omega); %arcseconds

%mean obliquity of the ecliptic
epsilon_del_arcsec = -46.8150*T - 0.00059*(T.^2) + 0.001813*(T.^3); %(Equation 22.2 in Meeus (1998))
epsilon = 23 + 26/60 + 21.448/3600 + epsilon_del_arcsec/3600;

%nutation in obliquity
delta_epsilon_arcsec = 9.20*cosd(omega) + 0.57*cosd(2*L) + 0.10*cosd(2*L_prime) - 0.09*cosd(2*omega);
%true obliquity, corrected for nutation according to Ch. 22, pg. 147
epsilon = epsilon + delta_epsilon_arcsec/3600; 

delta_So = delta_psi_arcsec.*cosd(epsilon); %in arcseconds
So = So + delta_So/3600;

%Finally add the longitude of the observer to the sidereal time at Greenwich
s = So + long; % (pg 92 in Ch. 13, sign switched because I define longitude negative in the Western Hemispherte)

So = mod(So,360); 
So = So/15; 

s = mod(s,360); 
s = s/15; 