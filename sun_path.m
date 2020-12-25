function sun_path(lat, dec)
%Function sun_path: calculates and plots the path of the Sun on the
%Celestial Sphere, given latitude of the observer and solar declination.
%
%Inputs:
%lat - latitude of the observer
%dec - declination of the Sun
%  Either lat or dec can be vectors, but not both! Function will plot
%  multiple solar paths in that case, for each element of the vector
%  variable.
%  It is not recommended for lat or dec to have more than ~10 elements, as
%  this will overwhelm the plot.
%  Allows for solar declination to exceed Earth's J2000 obliquity for
%  demo/pedagogical reasons; gives a warning if it does
%
%Author: Dr. Tihomir S. Kostadinov, 2014-2020 (major improvements in 2019-2020, e.g. addition of multiple path plotting)

obliquity = 23 + 26/60 + 21.448/3600; %Eq. 22.2 in Meeus (1998) for J2000.0, in degrees
if nargin==0
    %Provide default values
    lat = 43;
    dec = obliquity;
end

lat = lat(:);
dec=dec(:);

%Check that lat & dec are within allowed values & not both vectors
if ~(isnumeric(lat) && isnumeric(dec) && ~isempty(lat) && ~isempty(dec) && all(abs(lat)<=90) && all(abs(dec)<=90) )
    error('Incorrect input')
end
if (numel(lat)>1 && ~isscalar(dec)) || (numel(dec)>1 && ~isscalar(lat))
    error('Either lat or dec can be vectors, but not both!')
end

if any(abs(dec)>obliquity)
    warning('Some declination values exceed the J2000 obliquity of the Ecliptic in absolute value. Sun paths will be plotted, but they do not represent possible paths on contemporary Earth.');
end

%Convert inputs to radians
lat = lat*(pi/180);
dec = dec*(pi/180);

%Prepare figure for plotting
f1h = figure(1);
clf;

set(f1h,'Name','Daily path of the Sun on Celestial Sphere','NumberTitle','off','Units','normalized','OuterPosition',[0.3 0.15 .7 .85]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   PLOT THE CELESTIAL SPHERE     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 1;
[phi,theta] = meshgrid([0:10:180, 180]*pi/180,[0:10:360,360]*pi/180);

%parameterize the sphere for 3D plotting;
x = R*sin(phi).*cos(theta);
y = R*sin(phi).*sin(theta);
z = R*cos(phi);

h = mesh(x,y,z,ones(size(z)));
caxis([0 1.2])
colormap('gray')
alpha(0.3);
set(h,'EdgeColor','interp');
set(h,'FaceAlpha',2/3)
set(gca, 'DataAspectRatio', [ 1 1 1 ]);
hold on
shading interp
view(160,25)

%Hour angle of the true Sun, t
t = [0:.01:2*pi 2*pi]'; %roughly corresponds to the passage of time
%t = [0:.01:pi/2]'; %Allows testing of progress of parameterization from HA = 0

%plot the horizon:
horizon = plot3(sin(t),cos(t),zeros(size(t)),'b-','LineWidth',1.5);

%Plot the local meridian
meridian = plot3(sin(t),zeros(size(t)),cos(t),'g-','LineWidth',1.5);

%plot the local zenith axis:
h = line([0 0],[0 0],[-1.3 1.3]);
set(h,'Color','k')
text(0, 0, 1.4,'Zenith','FontSize',14);
text(0, 0, -1.4,'Nadir','FontSize',14);

%plot the local N-S axis:
h = line([-1 1],[0 0],[0 0]);
set(h,'Color','k')
text(1.05, 0, 0,'S','FontSize',14);
text(-1.05, 0, 0,'N','FontSize',14);

%plot the local E-W axis:
h = line([0 0],[-1 1],[0 0]);
set(h,'Color','k')
text(0, 1.05,0,'E','FontSize',14);
text(0, -1.05,0,'W','FontSize',14);

leg{1} = 'Local Horizon';
leg{2} = 'Local Meridian';

%Calculate alt/az solar coordinates as a function of hour angle
%then calculate Cartesian coordinates of Sun

%Case 1: plot for many declinations (or if both lat & dec are scalars):
if isscalar(lat)
    cm = autumn(numel(dec));
    
    %plot the polar axis towards the North & South Celestial Poles; This line coincides with the axis of rotation of Earth
    poles = line(1.15*[sin(pi/2-lat) -sin(pi/2-lat)],[0 0],1.15*[-cos(pi/2-lat) cos(pi/2-lat)],'LineWidth',1.5,'LineStyle','--');
    set(poles,'Color','c')
    text(-1.2*sin(pi/2-lat),0, 1.2*cos(pi/2-lat),'N pole','FontSize',14);
    text(1.2*sin(pi/2-lat),0,-1.2*cos(pi/2-lat),'S pole','FontSize',14);
    
    sun_path_h = NaN(numel(dec),1);
    for i =1:numel(dec)
        h = asin(sin(lat)*sin(dec(i)) + cos(lat)*cos(dec(i))*cos(t)); %Meeus (1998), Eq. 13.6
        A = atan2(sin(t),cos(t)*sin(lat) - tan(dec(i))*cos(lat)); %Meeus (1998), Eq. 13.5 -
        %AZIMUTH is this way measured from the SOUTH, which is more convenient for plotting below!
        
        %covert to Cartesian for plotting:
        % Astronomical (h,A) differ from the mathemaical spherical coordinates (theta,phi), so transfrom between the
        % two first (Phi is measured from the north pole (+ve z axis) in
        % mathematics, h is measured from the Equator in astronomy)
        % Azimuth is measured CW in astronomy, corresponding theta angle is measured CCW in mathematics
        x = R*sin(pi/2-h).*cos(-A); %This ensures correct progress of the Sun path with hour angle, positive westwards from local meridian
        %%%%%
        %Idea: consider animating the progress of the Sun on its path!
        %%%%%
        y = R*sin(pi/2-h).*sin(-A);
        z = R*cos(pi/2-h);
        
        sun_path_h(i) = plot3(x,y,z,'Color',cm(i,:),'LineWidth',1.5);
        leg{i+3} = ['Sun path, solar \delta = ',sprintf('%3.2f',dec(i)*180/pi),'^o'];
    end
    leg{3} = 'Earth''s axis of rotation ';
    title(['Daily path of the Sun in the sky at a latitude of ',sprintf('%3.2f',lat*180/pi),'^o'],'FontSize',14);
else
    %Case 2: plot for many latitudes:
    cm = autumn(numel(lat));
    
    %plot the polar axis towards the North & South Celestial Poles; This line coincides with the axis of rotation of Earth
    set(gca,'ColorOrder',autumn(numel(lat)));
    poles = line(1.15*[sin(pi/2-lat) -sin(pi/2-lat)],[0 0],1.15*[-cos(pi/2-lat) cos(pi/2-lat)],'LineWidth',1.5,'LineStyle','--');
    text(-1.2*sin(pi/2-lat),zeros(size(lat)), 1.2*cos(pi/2-lat),'N pole','FontSize',10);
    text(1.2*sin(pi/2-lat),zeros(size(lat)),-1.2*cos(pi/2-lat),'S pole','FontSize',10);
    
    sun_path_h = NaN(numel(lat),1);
    for i =1:numel(lat)
        h = asin(sin(lat(i))*sin(dec) + cos(lat(i))*cos(dec)*cos(t)); %Meeus (1998), Eq. 13.6
        A = atan2(sin(t),cos(t)*sin(lat(i)) - tan(dec)*cos(lat(i))); %Meeus (1998), Eq. 13.5 -
        
        x = R*sin(pi/2-h).*cos(-A);
        y = R*sin(pi/2-h).*sin(-A);
        z = R*cos(pi/2-h);
        
        sun_path_h(i) = plot3(x,y,z,'Color',cm(i,:),'LineWidth',1.5);
        leg{i+3} = ['Sun path, latitude \phi = ',sprintf('%3.2f',lat(i)*180/pi),'^o'];
    end
    leg{3} = 'Earth''s axis of rotation (one per latitude)';
    title(['Daily path of the Sun in the sky for solar declination \delta = ', sprintf('%3.2f',dec*180/pi),'^o'],'FontSize',14);
end

xlabel('N-S direction')
ylabel('E-W direction')
zlabel('Up-Down direction')
legend([horizon meridian poles(1) sun_path_h'],leg,'FontSize',12);
axis off