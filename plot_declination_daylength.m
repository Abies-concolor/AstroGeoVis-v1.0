function [MM,DD,DofY,Sun_delta_mu, Sun_delta_stdev, Sun_delta_rng, daylength] = plot_declination_daylength()
%Calculate and plot solar declination and duration of daylight for a generic representative year of 365 days. 
%Using actual declination from 4 consecutive years. 
%
%OUTPUT: 
%  MM - month
%  DD - day
%  DofY - day of year
%  Sun_delta_mu - average solar declination across 4 representative years, for each DofY
%  Sun_delta_stdev - average solar declination across these 4 representative years, for each DofY
%  Sun_delta_rng  - range of solar declination across these 4 representative years, for each DofY
%  daylength - daylength in hours on a lat/time grid given by meshgrid([1:365],[-90:0.1:90]);

%Author: Dr. Tihomir Kostadinov, 2018-2020
%Declination averagring as in Kostadiov and Gilb (2014), GMD, i.e. removing the leap day

t = datetime(2020,1,1,12,0,0,0):1:datetime(2023,12,31,12,0,0,0);
t(t==datetime(2020,2,29,12,0,0,0)) =[]; %Remove Feb. 29th from 2020

MMDD = dofyear2date_v2020(2021,1:365);
MM = MMDD(:,1); DD = MMDD(:,2);
DofY = (1:365)';

JD = date2jd_vec(t.Year,t.Month,t.Day,t.Hour,t.Minute,t.Second,'G');
[~, Sun_delta, ~] = solar_coord(JD);

Sun_delta = reshape(Sun_delta,365,4);
Sun_delta_mu = mean(Sun_delta,2);
Sun_delta_stdev = std(Sun_delta,[],2);
Sun_delta_rng = max(Sun_delta,2)-min(Sun_delta,2);

%%% Daylength calculations on a lat-time grid:
[xx,yy] = meshgrid([1:365],[-90:0.1:90]); %Grid in days of year and degrees of latitude
daylength = real(2/15*acosd(-tand(Sun_delta_mu(xx)).*tand(yy))); %In hours

%Daylength & declination figure for declinations approximating a 'generic' year of 365 days:
figure('Name','Global Daylength Patterns','NumberTitle','off','Units','normalized','OuterPosition',[.1 .1  .85 .85]);
yyaxis left
imagesc(xx(1,:),yy(:,1),daylength)
axis xy
axis([1 365 -90 90])
caxis([0 24])
pp = colorbar;
pp.Ticks = [0:24]; 
pp.Label.String = 'Duration of daylight, hours';
set(gca,'FontSize',18)
set(pp,'FontSize',14)
title('Duration of Daylight & Latitude of Subsolar Point')
ylabel('Latitude, degrees')
xlabel('Date')
set(gca,'YTick', [-90:15:90])
set(gca,'XTick',dofyear2date_v2020(2021,[1:12]',[15 15 20 15 15 21 15 15 22 15 15 21]')) %any non-leap year will do
set(gca,'XTickLabel',{'Jan 15','Feb 15', 'Mar 20 \newlineEquinox', 'Apr 15', 'May 15', 'Jun 21 \newlineSolstice', 'Jul 15', 'Aug 15',...
    'Sep 22 \newlineEquinox', 'Oct 15', 'Nov 15', 'Dec 21 \newlineSolstice'})
set(gca,'FontSize',14)
hold on
plot([1 365],[0 0],'k-.','LineWidth',1.2);
plot([dofyear2date_v2020(2021,3,20) dofyear2date_v2020(2021,3,20)],[-90 90],'-','LineWidth',1.2,'color',[.85 .85 .85]);
plot([dofyear2date_v2020(2021,6,21) dofyear2date_v2020(2021,6,21)],[-90 90],'-.','LineWidth',1.2,'color',[.85 .85 .85]);
plot([dofyear2date_v2020(2021,9,22) dofyear2date_v2020(2021,9,22)],[-90 90],'-','LineWidth',1.2,'color',[.85 .85 .85]);
plot([dofyear2date_v2020(2021,12,21) dofyear2date_v2020(2021,12,21)],[-90 90],'-.','LineWidth',1.2,'color',[.85 .85 .85]);

yyaxis right
plot(DofY,Sun_delta_mu,'LineWidth',1.75) %,)
ylabel('Latitude of Subsolar Point, degrees')
%Students can investigate the rate of change of declination and add it to
%the plot. Note discontinuity due to the choice of how leap years were treated. 
end