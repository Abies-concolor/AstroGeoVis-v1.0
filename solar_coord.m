function [RA_Sun, dec_Sun, Rvector] = solar_coord(JD)
%Outputs the solar RA and declination, and the Earth's radius-vector length, 
%at the given instant of JD, according to
%Meeus (1998), Astronomical Algorithms, Ch. 25 - the Low accuracy method
%Implemented & vectorized for MATLAB(r) by Dr. T. S. Kostadinov

%INPUT
% JD - Julian day number, provide double number for higher accuracy; 
%      (refers to an instant of time referenced to Universal Time(UT))
%      See date2jd_vec.m and Meeus(1999) for more details. 
%      Nx1 vector, where N is the number of observations. 

%OUTPUT: 
% RA_Sun - Apparent right ascension of the Sun, in degrees
% dec_Sun - Apparent declination of the Sun, in degrees
% Rvector - Sun-Earth radius vector length, in Astronomical units (AU)

if ~strcmpi(class(JD),'double')
    warning('Input JD should be a double for sufficient accuracy to be achieved');
end

JD = JD(:);
%The delta-T correction in the calculation of JD here has been ignored, see Meeus (1998), Ch. 10
T = (JD - 2451545.0)./36525;

Lo = 280.46646 +36000.76983*T + 0.0003032*(T.^2); %geometric mean longitude of the Sun, referred to the mean equinox of date
M = 357.52911 + 35999.05029*T - 0.0001537*(T.^2); %Mean anomaly of the Sun (same as mean anomaly of Earth)(Eq. 25.3)
%Note that a slightly different expression for M is given in Ch. 22, pg. 144, here I retain the expression in Ch. 25 given for low accuracy 
%The difference between the two formulas does not exceed 5.4 arcsec over 2
%centuries before and after J2000, smaller than the accuracy of 0.01 degrees for Ch. 25 stated for the low accuracy formulas

e = 0.016708634 - 0.000042037*T - 0.0000001267*(T.^2); %eccentricity of Earth's orbit

C = (1.914602 - 0.004817*T - 0.000014*(T.^2)).*sind(M) + (0.019993 - 0.000101*T).*sind(2*M) + 0.000289.*sind(3*M);

lambda_true = Lo + C; %Meeus Sun symbol, true longitude of the Sun
nu = M + C; %Meeus nu, Sun true anomaly

Rvector = (1.000001018)*(1-e.^2)./(1+e.*cosd(nu));

%apparent longitude lambda, corrected for nutation and aberration
%omega = 125.04 - 1934.136.*T;
omega = 125.04452 - 1934.136261.*T; %Omega expression in Ch. 22, pg 144 of Meeus (1998)
lambda = lambda_true - 0.00569 - 0.00478*sind(omega);

%mean obliquity of the ecliptic
epsilon_del_arcsec = -46.8150*T - 0.00059*(T.^2) + 0.001813*(T.^3); %(Equation 22.2 in Meeus (1998))
epsilon = 23 + 26/60 + 21.448/3600 + epsilon_del_arcsec/3600;

%Calculation of nutation in obliquity according to Ch. 22 (pg. 144), approximated by 0.00256*cosd(omega) in Ch. 25 (low accuracy);
L_prime = 218.3165 + 481267.8813*T; %Mean longitude of the Moon
delta_epsilon_arcsec = 9.20*cosd(omega) + 0.57*cosd(2*Lo) + 0.10*cosd(2*L_prime) - 0.09*cosd(2*omega);
epsilon = epsilon + delta_epsilon_arcsec/3600; 

%apparent RA and dec of Sun
RA_Sun = atan2d(cosd(epsilon).*sind(lambda),cosd(lambda));
qq = RA_Sun <0; 
RA_Sun(qq) = RA_Sun(qq) + 360;

dec_Sun = asind(sind(epsilon).*sind(lambda));