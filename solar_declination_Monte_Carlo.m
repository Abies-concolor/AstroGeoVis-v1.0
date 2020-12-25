%Monte Carlo simulation for analysis and propagation of error for the solar declination exercise:
%Author: Dr. Tihomir S. Kostadinov, 2012-2020

%Note that analytical propagation of error via derivatives can also be used.
 
%In the Monte Carlo analysis below, we consider the other inputs constants, known accurately with very
%high precision - e.g. the latitude of observer, PVC pipe length. 

%Do a simple Monte Carlo simulation to propagate error in one measurement -
%say shadow length, to the final result - declination of the Sun. Assume a
%normal distribution of errors around a certain true value, see what
%happens to simulated declination calculations - ask students
%to produce histogram of errors & their stats. --> for more advanced coding, math/(geospatial)stats or (geo)physics students. 

%Let the exact azimuth be 250 degrees and the exact shadow length 175 cm.
%Now assume error with normal distribution of mean zero and stdev = 3
%degrees for the azimuth measurement; and error with normal distribution
%with zero mean and stdev of 1.5 cm for the shadow length measurement. 
%Simulate 1,000,000 such measurements and show a histogram of the final
%declination computations.  How has the error in measurements propagated to
%declination? 

lat = 33; %deg
pipe_length = 61; %cm

%First vary the shadow length
A = 250;%degrees
shadow_length = 175 + randn(10^6,1)*1.5; %centimeters

h = atand(pipe_length./shadow_length);
delta = asind(sind(h)*sind(lat) + (cosd(h)*cosd(lat)).*cosd(A)); 

h_actual = atand(pipe_length/175);
delta_actual = asind(sind(h_actual)*sind(lat) + (cosd(h_actual)*cosd(lat)).*cosd(A)); 

mean(delta)
std(delta)
min(delta)
max(delta)
median(delta)

figure
histogram(delta)
alim = axis; 
hold on
plot([delta_actual delta_actual],[alim(3) alim(4)],'r-','LineWidth',2)
xlabel('Computed solar \delta, deg')
ylabel('Histogram bin count')
title('Shadow length error propagation to solar \delta')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then vary the azimuth
A = 250+ randn(10^6,1)*3;%degrees
shadow_length = 175;  %centimeters

h = atand(pipe_length./shadow_length);
delta = asind(sind(h)*sind(lat) + (cosd(h)*cosd(lat)).*cosd(A)); 

A_actual = 250;
delta_actual = asind(sind(h)*sind(lat) + (cosd(h_actual)*cosd(lat)).*cosd(A_actual)); 

mean(delta)
std(delta)
min(delta)
max(delta)
median(delta)

figure
hh = histogram(delta)
alim = axis; 
hold on
plot([delta_actual delta_actual],[alim(3) alim(4)],'r-','LineWidth',2)
xlabel('Computed solar \delta, deg')
ylabel('Histogram bin count')
hh.EdgeColor = [0.1 0.4 0.9];
title('Azimuth error propagation to solar \delta')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%In this case it turns out the result is NOT very sensitive to the
%shadow length, but it IS quite sensitive to the azimuth measurement. Whether this is 
%the case in general can be ascertained with derivative analysis; of course
%this also depends on the assumed input variable error statistics. 