function [r,theta] = gnomon_trace(lat, dec)
%Calculates and plots the trace of the shadow tips, cast on a horizontal
%surface by a vertical pole (GNOMON, stylus) of length 1 m, positioned at the plot origin.
%
%INPUTS: lat, Latitude of the observer and dec, solar declination, both in degrees 
%OUTPUTS: A polar plot as described above and the polar coordinates (r and theta pairs) used to produce the plot
%
%Plot may not work well at high latitudes
%Author: Dr. Tihomir S. Kostadinov - 2014-2020

obliquity = 23 + 26/60 + 21.448/3600; %Eq. 22.2 in Meeus (1998) for J2000.0

if nargin==0
    %Provide default values 
    lat = 33+7/60+45/3600;
    dec = obliquity;
end

%Convert inputs to radians
lat = lat*(pi/180); 
dec = dec*(pi/180); 

%Hour angle of the Sun, denoted t here, is used as the time parameter - time steps in true Sun hour angle are 
%roughly equivalent to time steps in civilian time over the course of a
%day, the small difefernces are due to the Equation of Time. See analemma.m
t = [-pi:0.01:pi, pi]';
h = asin(sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(t)); %Meeus (1998) Eq. 13.6
A = atan2(sin(t),cos(t)*sin(lat) - tan(dec)*cos(lat)) + pi; %Meeus (1998) Eq.13.5

r = cot(h);
theta = (pi/2)-A + pi; 
qq = theta<0; 
theta(qq) = theta(qq) + 2*pi;

%remove all values where h <= 12 degrees - shadows are too long and plot doesn't look good
q = h<=12*pi/180; 
r(q) = [];
theta(q) = [];
t_truncated = t;
t_truncated(q) = [];

HAp12 = t_truncated*(180/pi)/15 + 12;

%Polar plot
fh = figure; 
set(fh,'Name','Daily path of the shadow of the Sun on a level surface, cast by 1 m vertical pole','NumberTitle','off','Units','normalized','OuterPosition',[0.3 0.15 .7 .85]);
polarscatter(theta,r,65,HAp12,'.')
thetaticknums = 0:30:330;
thetaticks(thetaticknums);  
p = {'120^o','150^o','180^o (South)','210^o','240^o','270^o (West)','300^o','330^o','0^o (North)','30^o','60^o','90^o (East)'};
thetaticklabels(flipud(p'))
ww = colorbar;
mver = ver('MATLAB');
if strncmp(mver.Version,'9.9',3) %new turbo colormap will not work in recent older versions
   colormap turbo
else
   colormap jet
end
ww.Ticks = ceil(min(HAp12)):floor(max(HAp12));
ww.Label.String = 'Hour Angle of the Sun + 12 hrs';
ww.Label.FontSize = 12;
title(['Vertical Gnomon Shadow Tip Trace for \phi = ',num2str(lat*180/pi),'^o and \delta = ',num2str(dec*180/pi),'^o' ])
set(gca,'FontSize',14)

%Question to explore - why is the shadow path always a straight line for a declination
%of exactly 0 degrees? - Solve for this analytically