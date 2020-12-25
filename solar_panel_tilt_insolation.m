function [Energy_produced_kWh,panel_insolation] = solar_panel_tilt_insolation(lat,lon,YYYY,MM,DD,TZ,alpha_p,A_p,Pmax,dt,to_plot)
% Solar Panel Tilt and Insolation Analysis for a single day
%
%%% INPUTS:
%% lat - Latitude of solar panel location, degrees
%% lon - Longitude of the solar panel location, degrees (negative in Western Hemisphere)
%% YYYY - year, if empty, 2020 is assumed, and the Gregorian calendar is assumed
%% MM - month
%% DD - day
%% TZ - time zone offset in hours, from UT, negative in Western Hemisphere
%%%%%%% Solar panel position inputs:
%% alpha_p - Tilt of the solar panel, in degrees. Defined as angle between the normal to the plane of the solar panel and the local vertical direction.
%% A_p - Azimuth of the panel orientation, in degrees, clockwie from true North. This is the azimuth of the projection of the normal to the plane of the
%%       solar panel on the plane of the local horizon. May result in degeneracy for case of panel facing towards zenith.
%% Pmax - the maximum power at STC for the solar panel, typically given in panel specs. In Watts.
%% dt - time step for calculations in seconds
%% to_plot - flag whether to plot time series of daily solar & panel geometry figures. 1 - plot, 0 - do not plot.
%%           set to 0 if this function is called multiple times in a loop.
% Author: Dr. Tihomir Kostadinov, CSUSM, Sept. 2019 - Dec. 2020

%Default values - location is CSU San Marcos campus:
%Note - rigorous input checking is not done.
if nargin==0
    lat = 33+7/60+45/3600;
    lon = -117 -9/60 - 30/3600;
    YYYY = 2020;
    MM = 10;
    DD = 25;
    TZ = -8;
    alpha_p = lat; %Panel tilted at latitude, i.e. panel is perpendicular to the plane of the Equator
    A_p = 180; %the default for extratropical North Hemisphere deployments, panel facing southward
    Pmax = 100;
    dt = 60;
    to_plot = 1;
end
%Allow not specifying some arguments - input as empty to use these defaults:
if isempty(lat)
    lat = 33+7/60+45/3600;
end
if isempty(lon)
    lon = -117 -9/60 - 30/3600;
end
if isempty(YYYY)
    YYYY = 2020;
end
if isempty(MM)
    MM = 10;
end
if isempty(DD)
    DD = 25;
end
if isempty(TZ)
    TZ = -8;
end
if isempty(alpha_p)
    alpha_p = lat; %Panel tilted at latitude, i.e. panel is perpendicular to the plane of the Equator
end
if isempty(A_p)
    A_p = 180; %the default for extratropical North Hemisphere deployments, panel facing southward
end
if isempty(Pmax)
    Pmax = 100;
end
if isempty(dt)
    dt = 60;
end
if isempty(to_plot)
    to_plot = 1;
end

%Constants
So = 1361; %W/m^2, TOA irradiance at 1 AU and perpendicular to solar rays - See Kopp and Lean (2011), GRL, doi:10.1029/2010GL045777
standard_irradiance = 1000; %W/m^2; The irradiance for which Pmax is given.
%Typical value used in STC, e.g. from panel specs of the (c) Renogy RNG-100D 100-Watt panel

%Assuming local standard time, ofset by TZ from UT, TZ given above
%Calculate coordinates of the Sun for small time steps for the entire 24 hr
%period. Periods when the Sun is below the horizon or not shining on panel
%will be excluded from the calculation later
t = datetime(YYYY,MM,DD,0,0,0):dt/86400:datetime(YYYY,MM,DD,23,59,59); %time in steps of dt seconds in local civilian time
t = t';
t_local = t;
t = t-TZ/24; %t is now in UT

%Julian Day
JD = date2jd_vec(t.Year,t.Month,t.Day,t.Hour,t.Minute,t.Second,'G');%Time has to be Universal time!!!
[RA_Sun, delta_Sun, Rvec] = solar_coord(JD);
[s,~] = sidereal_time(JD,lon);

%Convert to local horizontal coordinates
H = s - RA_Sun/15; %H is in hours here, s in hours, RA_Sun in degrees
h_Sun = asind(sind(delta_Sun)*sind(lat)+ cosd(delta_Sun)*cosd(lat).*cosd(H*15) );
A_Sun = atan2d(sind(H*15),cosd(H*15)*sind(lat) - tand(delta_Sun)*cosd(lat)); %Meeus 1998, Eq. 13.5, page 93
A_Sun = A_Sun + 180; %convert azimuth to degrees clockwise from North

%Plot local horizontal coordinates of the Sun, investigate maximum altitude vs. local solar noon timing
%Find and indicate local solar noon - transit time, compare with timing of maximum altitude
[h_max,idx] = max(h_Sun);
[~,idx2] = min(abs(A_Sun-180));
[H_min,idx3] = min(abs(H));

%Note that
idx2==idx3;
%indicating consistency, i.e. local transit occurs when H = 0 h, and A = 180 deg.,
%however, the highest altitude h_max occurs 11 seconds earlier
%Maximum altitude:
t_local(idx);
t_local(idx2);
%This occurs because declination is also changing continuously in the Meeus equations used here. For a fixed declination maximum altitude will occur at trasnit as expected.

%Calculate angle between direction of the Sun and the normal to the solar panel plane
%Direction of the Sun as unit Cartesian vector in (h,A) coordinate system
Sun_dir_XYZ = [cosd(h_Sun).*sind(A_Sun), cosd(h_Sun).*cosd(A_Sun), sind(h_Sun)];
%Direction of the normal to the solar panel as unit Cartesian vector in (h,A) coordinate system
panel_normal_dir_XYZ = [sind(alpha_p)*sind(A_p), sind(alpha_p)*cosd(A_p), cosd(alpha_p) ];
%Dot product of these two vectors will give cosine of angle between them, i.e. sza_p, the solar zenith angle as measured in the plane of the solar panel
cos_sza_p = dot(Sun_dir_XYZ,repmat(panel_normal_dir_XYZ,size(Sun_dir_XYZ,1),1),2);
sza_p = acosd(cos_sza_p); %This is simplified because the vectors involved have norm 1 (i.e. are unit vectors).

%Insolation at the solar panel - corrected for Sun-Earth distance is thus:
panel_insolation = So*cosd(sza_p).*((1./Rvec).^2);
%No atmospheric effects, TOA irradiance assumed, and no atmospehric
%refraction correction is implemented. No local shading effects taken into
%account.

horizontal_insolation = So*sind(h_Sun).*((1./Rvec).^2); %for comparison
horizontal_insolation(horizontal_insolation<0) = 0; %Sun below horizon

%Sun not shining on panel surface (below the horizon or shining on the rear surface of the panel) - indicated by negative
%cos_sza_p values and/or negative h_Sun values, these cases should be removed and not used, replaced by zeros - they mean insolation 0 W/m^2
qq = panel_insolation<0 | h_Sun <= 0;
panel_insolation(qq) = 0;

%Scale panel wattage with actual irradiance, assume otherwise conditions similar to standard
panel_wattage = Pmax*(panel_insolation/standard_irradiance); %For simplicity assumes Pmax can be exceeded if irradiance is higher than standard

Energy_produced_J = trapz(panel_wattage)*dt; %multiply by time step, so this integration yields Joules
Energy_produced_kWh = Energy_produced_J/3.6e+6;

if to_plot
    figure
    yyaxis left
    plot(t_local,h_Sun,'.')
    hold on
    plot([t_local(1) t_local(end)],[0 0],'k--')
    ylabel('Altitude of the Sun, deg')
    yyaxis right
    plot(t+TZ/24,A_Sun,'.')
    ylim([0 360])
    plot([t(idx2) t(idx2)]+TZ/24,[0 360],'b--')
    xlabel('Date and time, PST (UT-8)')
    ylabel('Azimuth of the Sun, deg CW from true North')
    set(gca, 'XTick', datetime(YYYY,MM,DD,0,0,0):3/24:datetime(YYYY,MM,DD,24,0,0))
    title(['Horizontal Coordinates of the Sun for ', datestr(t_local(1),1), '; \phi = ', num2str(lat),...
        '^o, \lambda = ', num2str(lon),'^o'])
    
    figure
    yyaxis right
    plot(t_local, panel_insolation)
    hold on
    plot(t_local, horizontal_insolation,'--')
    ylabel('Irradiance, W m^-^2')
    yyaxis left
    plot(t_local,sza_p)
    ylim([0 180])
    set(gca,'YTick',[0:45:180])
    xlabel('Date and time, PST (UT-8)')
    ylabel('Angle between direction to Sun and normal to panel plane, deg')
    legend({'Panel-Sun angle, deg','Panel irradiance, W m^-^2','Horizontal Irradiance, W m^-^2',})
    set(gca, 'XTick', datetime(YYYY,MM,DD,0,0,0):3/24:datetime(YYYY,MM,DD,24,0,0))
    title(['Panel-Sun Geometry & Irradiances on ', datestr(t_local(1),1), ' at \phi = ', num2str(lat),...
        '^o, \lambda = ', num2str(lon),'^o \newline                                          ', ...
        'Total energy produced: ',num2str(Energy_produced_kWh), ' kWh']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3D plotting of the Sun and panel orientation geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TK, September 29, 2019 - Oct. 27, 2020
    
    % Normal vector to the plane of the solar panel - from above - scale for visibility:
    N = 1.15*panel_normal_dir_XYZ';
    
    figure
    %The plane passing thorugh the origin with (A,B,C) as a normal vector has
    %equation Ax = By + Cz = 0, see e.g. https://www.math.tamu.edu/~glahodny/Math251/Section%2011.4.pdf
    [x,y] = meshgrid([-0.5:0.01:0.5],[-0.5:0.01:0.5]);
    z = -(N(1)/N(3))*x - (N(2)/N(3))*y;
    mesh(x,y,z)
    hold on
    
    quiver3(0,0,0,N(1),N(2),N(3),0,'LineWidth',2)
    set(gca,'DataAspectRatio', [ 1 1 1]);
    xlabel('East-West'), ylabel('North-South'),zlabel('Up-Down')
    
    %Now plot several vectors through the origin that are the direction of the Sun at various times, from sunrise to sunset only:
    pp = find(h_Sun >=0);
    if ~isempty(pp)
        Svec_to_plot = Sun_dir_XYZ([pp(1:3600/dt:end); pp(end)],:); %Sunrise to sunset in 1 hr steps
        Svec_LSN = 1.3*Sun_dir_XYZ(idx3,:); %Local solar noon
        quiver3(zeros(size(Svec_to_plot,1),1),zeros(size(Svec_to_plot,1),1),zeros(size(Svec_to_plot,1),1),...
            Svec_to_plot(:,1),Svec_to_plot(:,2),Svec_to_plot(:,3),0,'color',[1 0.8 0.1],'LineWidth',2)
        quiver3(0,0,0,Svec_LSN(1),Svec_LSN(2),Svec_LSN(3),0,'color',[.75 0.1 0.1],'LineWidth',2);
        legend({'Plane of solar panel','Normal to solar panel plane','Direction to Sun','Direction to Sun @transit'})
    else
        warning('The Sun is never above the horizon!')
        legend({'Plane of solar panel','Normal to solar panel plane'})
    end   
    view(-60,10) %may not be suitable for all gemoetries, suitable for default values.
    %End 3D plotting of Sun-panel geometry
end