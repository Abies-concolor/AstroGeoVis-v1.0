function [Energy_annual_TS, Energy_annual_total] = solar_panel_annual(lat,lon,YYYY,TZ,alpha_p,A_p,Pmax,dt,to_plot)
%%% INPUTS:
%% lat - Latitude of solar panel location, degrees
%% lon - Longitude of the solar panel location, degrees (negative in Western Hemisphere)
%% YYYY - year, if empty, 2020 is assumed, and the Gregorian calendar is assumed
%% TZ - time zone offset in hours, from UT, negative in Western Hemisphere
%%%%%%% Solar panel position inputs:
%% alpha_p - Tilt of the solar panel, in degrees. Defined as angle between the normal to the plane of the solar panel and the local vertical direction.
%% A_p - Azimuth of the panel orientation, in degrees, clockwie from true North. This is the azimuth of the projection of the normal to the plane of the
%%       solar panel on the plane of the local horizon. May result in degeneracy for case of panel facing towards zenith.
%% Pmax - the maximum power at STC for the solar panel, typically given in panel specs. In Watts.
%% dt - time step for calculations in seconds (over a single day, used by solar_panel_tilt_insolation.m which is called here)
%% to_plot - flag whether to plot time series of daily solar panel energy production estimates; 1 - plot, 0 - do not plot.
%%
% Author: Dr. Tihomir Kostadinov, 2020

if nargin==0
    lat = 33+7/60+45/3600;
    lon = -117 -9/60 - 30/3600;
    YYYY = 2021;
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
    YYYY = 2021;
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

tt = datetime(YYYY,1,1,12,0,0):1:datetime(YYYY,12,31,12,0,0);

Energy_produced_kWh = NaN(size(tt));
for i = 1:numel(tt)
    [Energy_produced_kWh(i),~] = solar_panel_tilt_insolation(lat,lon,YYYY,tt(i).Month,tt(i).Day,TZ,alpha_p,A_p,Pmax,dt,0);
end

Energy_annual_TS = Energy_produced_kWh;
Energy_annual_total = sum(Energy_annual_TS);

if to_plot
    figure
    plot(tt,Energy_produced_kWh,'.-')
    hold on
    yl = ylim;
    %Approximate equinoxes and solstices, can be made exact by examining the
    %time series of solar declination produced in solar_panel_tilt_insolation.m
    plot([datetime(YYYY,3,20,12,0,0) datetime(YYYY,3,20,12,0,0)],[yl(1) yl(2)],'r-.','LineWidth',1.5);
    plot([datetime(YYYY,6,21,12,0,0) datetime(YYYY,6,21,12,0,0)],[yl(1) yl(2)],'r-.','LineWidth',1.5);
    plot([datetime(YYYY,9,22,12,0,0) datetime(YYYY,9,22,12,0,0)],[yl(1) yl(2)],'r-.','LineWidth',1.5);
    plot([datetime(YYYY,12,21,12,0,0) datetime(YYYY,12,21,12,0,0)],[yl(1) yl(2)],'r-.','LineWidth',1.5);
    set(gca,'XTick',datetime(YYYY,1:12,15,12,0,0))
    xtickformat('MMM.dd')
    xtickangle(45)
    xlabel('Date')
    ylabel('Total daily energy production, kWh')
    legend({'Daily panel energy production','Solstices/Equinoxes'})
    title(['Total annual energy production: ', sprintf('%5.2f',Energy_annual_total), '  kWh'])
end