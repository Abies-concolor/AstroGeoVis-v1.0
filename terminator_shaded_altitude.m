function terminator_shaded_altitude(Y,M,D,H,MN,S)
%Plot terminator on a map of the world in equidistant cylindircal projection
%Overlay a colormap of instantaneous solar altitude globally
%Display subsolar and anti-solar points.
%
%Author: Dr. T.S. Kostadinov, 2015-2020
%
%Inputs:
%Option 1: TZ: Enter your time zone offset from UTC, TZ, in hours, negative in the Western Hemisphere.
%Function will use your computer's local system date and time to plot the terminator
%
%Option 2: Enter exact date and time desired, in UT, Gregorian calendar assumed 
% Y - year
% M - month
% D - day
% H - hour (24 hr system)
% MN - minunte
% S - second; It is possible to skip entering the seconds, S, they are assumed equal to zero in that case (S= 0)

if nargin ==1
    %enter timezone offset from UTC, and assume user wants the current date/time to be used
    TZ = Y; %in hours, negative west of Greenwich
    dv = datevec(now - TZ/24); %Date/time in UTC
    Y = dv(1);
    M = dv(2);
    D = dv(3);
    H = dv(4);
    MN = dv(5);
    S = dv(6);
elseif nargin == 0
    %Example with UTC time given
    Y = 2020;
    M = 10;
    D = 7;
    H = 12;
    MN = 0;
    S = 0;
elseif nargin == 5
    S = 0;
elseif nargin ==6
    %do nothing, proceed
else
    error('Incorrect number of arguments')
end
%%%%%%%%%%%%%%%%%%%%

obliquity = 23 + 26/60 + 21.448/3600; %Eq. 22.2 in Meeus (1998) for J2000.0
%Value coincides within << 1 arcsec with the value used in Laskar et al.
%(2004) and adopted in the Earth orbit model of Kostadinov and Gilb (2014),
%GMD, namely 0.4090928042223415 radians, for J2000.0 
%Here obliquity is only used to plot the tropics on the map.

JD = date2jd_vec(Y,M,D,H,MN,S,'G'); %Julian day, input date/time has to be in UT!
[RA_Sun, dec_Sun, ~] = solar_coord(JD); %Apparent right ascension and declination of the Sun

[So, ~] = sidereal_time(JD,0); % Apparent sidereal time at the Prime Meridian

%Create a lat/long grid over which the gridded solar altitude computations will take place - 1 degree resolution:
[lon,lat] = meshgrid(-179:180,90:-1:-90);
s = NaN(size(lat,1),size(lat,2));

for i = 1:size(lat,1)
    for j = 1:size(lat,2)
        [s(i,j), ~] = sidereal_time(JD,lon(i,j)); %local sidereal time
    end
end

HA = s*15 - RA_Sun; %hour angle of the Sun in degrees
solar_alt = asind(sind(lat)*sind(dec_Sun)+cosd(lat)*cosd(dec_Sun).*cosd(HA));

%Ideas for additions/improvements:
%1) show a map of solar azimuths
%2) improve colormap
%3) show with distinct lines and/or color the different twilights
%4) solve for the time when a given azimuth of the Sun occurs throughout the year
%5) take into account atmospheric refraction and size of solar disk in altitude computations.

%Find lat/lon of subsolar point
subsolar_latitude = dec_Sun;
%So + longitude = local sidereal time = Sun right ascension when Sun is transiting (since hour angle = 0 hrs)
%Solve for longitude where that is true at the moment to find the subsolar
%longitude, i.e. the meridian where the Sun is transiting
subsolar_longitude = RA_Sun - 15*So; %needs to be in degrees

%Deal with cases when the subsolar longitude is not in -180 to 180 range - all cases should be treated exhaustively
if subsolar_longitude>180
    subsolar_longitude = subsolar_longitude - 360;
elseif subsolar_longitude<-180
    subsolar_longitude = subsolar_longitude + 360;
end
if subsolar_longitude > 0
    antisolar_longitude = subsolar_longitude - 180;
else
    antisolar_longitude = subsolar_longitude + 180;
end
antisolar_latitude = -subsolar_latitude;

% Parameterize the terminator as a great circle in a coordinate system in which it is the Equator with
% poles the subsolar and antisolar points, then change coordinate systsem to Earth's lat/lon.

%Parameterization of the Equator in Cartesian coordinates
t = [0:0.01:2*pi, 2*pi];
x = sin(t);
y = cos(t);
XX = [x', y', zeros(numel(t),1)];

%Tilt terminator so that it becomes the great circle with poles - the subsolar and the antisolar point, in Earth's lat/lon coordinate system
%This is using a rotation matrix about an axis. Script generate_rot_m.m is from Kostadinov and Gilb (2014), GMD.
tilt_m = generate_rot_m(deg2rad(90-subsolar_latitude), deg2rad(subsolar_longitude+90));
ZZ = tilt_m*XX';

%Convert the tilted terminator coordinates to lat/lon (modified spherical in essence)
phi = 90-acosd(ZZ(3,:));
theta = atan2d(ZZ(2,:),ZZ(1,:));

%For plotting purposes, rearrange theta and phi vectors so breakpoint between western and eastern hemispheres is at the beginning
% allows a line plot to be made without a weird line across hemispheres.
%The below may not work perfectly when the Sun has declination very close to zero.
%Alternate logic/approach should be explored.
qq = find(abs(diff(theta))>45);
if numel(qq)~=1
    warning('Breakpoint in longtiudes is not a single point along terminator points...')
    qq = qq(1);
end
theta = [theta(qq+1:end), theta(1:qq)];
phi = [phi(qq+1:end), phi(1:qq)];

%ready for plotting all elements
figure('Name','Terminator Plot & Solar Altitude Map','NumberTitle','off','Units','normalized','OuterPosition',[.1 .1  .85 .85]);
%plot solar altitude false color map
imagesc(lon(1,:),lat(:,1),solar_alt)
colormap(parula)
axis xy
hold on
caxis([-90 90])
colorbar
%Plot line of terminator - center of solar disk is at altitude 0 degrees along this line, without taking itno account atmospheric refraction
plot(theta,phi,'r.-','MarkerSize',6)
axis([-180 180 -90 90])

%Plot the coastlines, using the GSHHS coastlines data set.
try
    coastlines = load('GSHHS_2.3.7_hi_extracted_levels16.mat');
catch
    [b,a] = uigetfile('*.mat','Locate the coastlines data mat-file');
    coastlines = load(fullfile(a,b));
end

plot(coastlines.lon,coastlines.lat,'.','MarkerSize',.1,'Color',[.6 .6 .6])
grid on

%Plot subsolar and anti-solar points
plot(subsolar_longitude,subsolar_latitude,'r*','MarkerSize',24)
plot(antisolar_longitude,antisolar_latitude,'kx','MarkerSize',24)

%Plot meridians of true solar midnight and solar noon:
plot([subsolar_longitude, subsolar_longitude],[-90 90],'r-')
plot([antisolar_longitude, antisolar_longitude],[-90 90],'k-')

%Plot lines of tropics
plot([-180 180],[obliquity obliquity],'-.','Color',[.6 .6 .6])
plot([-180 180],[-obliquity -obliquity],'-.','Color',[.6 .6 .6])

set(gca,'XTick',[-150:30:150])
set(gca,'YTick',[-60:30:60])

mver = ver('MATLAB');
if strncmp(mver.Version,'9.9',3) %new subtitle command will not work in recent older versions
    title('Terminator (red line), subsolar (red star) and antisolar (black X) points, and solar altitude (deg) (color scale)');
    subtitle([datestr(datenum(Y,M,D,H,MN,S),31),' UT; \delta Sun = ',sprintf('%02.2f',dec_Sun),'^o' ]);
else
    title(['Terminator (red line), subsolar (red star) and antisolar (black X) points, and solar altitude (deg) (color scale) \newline',...
        '                                                    ',...
        datestr(datenum(Y,M,D,H,MN,S),31),' UT; \delta Sun = ',sprintf('%02.2f',dec_Sun),'^o' ]);
end
xlabel('Longitude, deg')
ylabel('Latitude, deg')
set(gca,'DataAspectRatio', [ 1 1 1 ]);