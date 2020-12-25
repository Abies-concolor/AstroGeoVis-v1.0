function [r,theta] = gnomon_trace_annual(lat,YYYY)
%Calculates and plots the trace of the shadow tips, cast on a horizontal
%surface by a vertical pole (GNOMON, stylus) of length 1 m, positioned at
%the plot origin, over the course of a year. 
%
%INPUTS: Latitude of observer, lat (degrees), and the desired year, YYYY
%OUTPUTS: A polar plot as described above and the r and theta pairs used to produce the plot
%
%Plot may not work well at high latitudes
%Author: Dr. Tihomir S. Kostadinov - 2014-2020

if nargin==0
    %Provide default values
    lat = 33+7/60+45/3600;
    YYYY = 2021;
end

long = -117 -9/60 - 30/3600;%positive East of Greenwich, degrees
TZ = -8; %positive East of Greenwich, hours offset from UT
longitude_correction = TZ - long/15; %in hours
fractional_hr = 12 - TZ + longitude_correction; %Time of local noon of the Mean Sun

t = datetime(YYYY,1,1,0,0,0):1:datetime(YYYY,12,31,0,0,0);
t = t';
t = t+fractional_hr/24;

JD = date2jd_vec(t.Year,t.Month,t.Day,t.Hour,t.Minute,t.Second,'G'); %Julian day
[RA_Sun, dec_Sun, ~] = solar_coord(JD); %RA & dec in degrees, Rvector in AU
[local_s, ~] = sidereal_time(JD,long); %sidereal time in hours;

H = local_s - RA_Sun/15; %(pg. 92 of Meeus (1998) - Ch. 13) - in hours
H = mod(H,24); %Add 24 to negative H values
H = H*15;

h = asin(sind(lat)*sind(dec_Sun) + cosd(lat)*cosd(dec_Sun).*cosd(H)); %Meeus(1998), Eq. 13.6
A = atan2(sind(H),cosd(H)*sind(lat) - tand(dec_Sun)*cosd(lat)) + pi; %Meeus(1998), Eq. 13.5

r = cot(h); %Stylus/gnomon length is considered unit (1 m), so it is not multiplied by here

theta = (pi/2)-A + pi;
qq = theta<0;
theta(qq) = theta(qq) + 2*pi;

%remove all values where h <= 0.1 degrees - shadows are too long and plot
%doesn't look good, also - negative h can occur near the poles at certain
%times of year and they would incorrectly plot with negative r's radially symmetrically
q = h<=0.1*pi/180;
r(q) = [];
theta(q) = [];

%Polar plot
fh = figure;
set(fh,'Name','Annual path of 1-m vertical gnomon shadow tip on a level surface','NumberTitle','off','Units','normalized','OuterPosition',[0.3 0.15 .7 .85]);
polarscatter(theta,r,65,[1:numel(r)],'.')
thetaticknums = 0:30:330;
thetaticks(thetaticknums);
p = {'120^o','150^o','180^o (South)','210^o','240^o','270^o (West)','300^o','330^o','0^o (North)','30^o','60^o','90^o (East)'};
thetaticklabels(flipud(p'))
%tticklabelnums = mod(90-thetaticknums,360); %Alternate way to convert from
%mathematical polar angle to geographical azimuth (but degree symbol harder to put)
ww = colorbar;
mver = ver('MATLAB');
if strncmp(mver.Version,'9.9',3) %new turbo colormap will not work in recent older versions
   colormap turbo
else
   colormap jet
end
ww.Ticks = dofyear2date_v2020(YYYY,1:12,ones(1,12));
ww.TickLabels = {'Jan. 1' 'Feb. 1' 'Mar. 1' 'Apr. 1' 'May 1' 'Jun. 1' 'Jul. 1' 'Aug. 1' 'Sep. 1' 'Oct. 1' 'Nov.1' 'Dec. 1'};
title(['Vertical Gnomon Shadow Tip Annual Trace for latitude = ',num2str(lat),'^o '])
set(gca,'FontSize',14)