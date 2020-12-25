%Plot global insolation for J2000 as calculated by the Earth Orbit v2.1
%model of Kostadinov and Gilb (2014), GMD, using Laskar et al. (2004) astronomical 
%solutions. Plot improved here for pedagogical purposes. See also plot_declination_daylength.m
%Author: Dr. Tihomir Kostadinov, Sept. 8, 2019 - Nov. 2020

insolationglobal = load('Global_Insolation_EOv2.1.dat');
%Data generated for So = 1361 W/m^2 and contemporary Laskar et al. (2004) solutions with Earth orbit model of Kostadinov and Gilb (2014) (GMD) 
So = 1361; %W/m^2 (Kopp and Lean (2011), GRL)
albedo = 0.3;
obliquity = 23 + 26/60 + 21.448/3600; %Eq. 22.2 in Meeus (1998) for J2000.0

QQ = insolationglobal(2:end,2:end);
dd = insolationglobal(1,2:end);
lat = insolationglobal(2:end,1);

figure('Name','Global Insolation','NumberTitle','off','Units','normalized','OuterPosition',[.1 .1  .85 .85]);
imagesc(dd,lat,QQ)
axis xy
title('Daily Average Solar Energy Receipt @ TOA, W m^-^2','FontSize',28)
xlabel('Date');
ylabel('Latitude,deg')
hold on
contour(dd,lat,QQ,[0:50:550],'k','ShowText','on')
set(gca,'XTick',dofyear2date_v2020(2017,[1:12]',[15 15 20 15 15 21 15 15 22 15 15 21]')) %any non-leap year will do
set(gca,'XTickLabel',{'Jan 15','Feb 15', 'Mar 20 \newlineEquinox', 'Apr 15', 'May 15', 'Jun 21 \newlineSolstice', 'Jul 15', 'Aug 15',...
    'Sep 22 \newlineEquinox', 'Oct 15', 'Nov 15', 'Dec 21 \newlineSolstice'})
set(gca,'FontSize',14)

plot([1 365],[0 0],'k-.','LineWidth',1.2);
plot([dofyear2date_v2020(2017,3,20) dofyear2date_v2020(2017,3,20)],[-90 90],'r-','LineWidth',1.2);
plot([dofyear2date_v2020(2017,6,21) dofyear2date_v2020(2017,6,21)],[-90 90],'r-.','LineWidth',1.2);
plot([dofyear2date_v2020(2017,9,22) dofyear2date_v2020(2017,9,22)],[-90 90],'r-','LineWidth',1.2);
plot([dofyear2date_v2020(2017,12,21) dofyear2date_v2020(2017,12,21)],[-90 90],'r-.','LineWidth',1.2);

plot([1 365],[obliquity obliquity],'-.','Color',[.3 .3 .3])
plot([1 365],[-obliquity -obliquity],'-.','Color',[.3 .3 .3])
plot([1 365],[90-obliquity 90-obliquity],'-.','Color',[.5 .5 .5])
plot([1 365],[obliquity-90 obliquity-90],'-.','Color',[.5 .5 .5])

ch = colorbar; 
ch.Label.String = 'Daily Insolation, W m^-^2';

%Mean annual insolation at each latitude:
pp = mean(QQ,2); %Mean across all days
pp_stdev = std(QQ,[],2);

figure('Name','Global Annual Insolation mean and standard deviation','NumberTitle','off','Units','normalized','OuterPosition',[.1 .1  .85 .85]);
yyaxis left
plot(lat,pp,'LineWidth',1.5)
axis([-90 90 135 450])
xlabel('Latitude, degrees')
ylabel('Average annual TOA insolation, W m^-^2')
title('Annual \mu and \sigma of TOA insolation by latitude, W m^-^2')
set(gca,'FontSize',18)
yyaxis right
plot(lat,pp_stdev,'x-','LineWidth',1.5)
axis([-90 90 0 250])
ylabel('Annual \sigma of TOA insolation, W m^-^2')

%Average insolation across all Earth (globally) and annually: 
%(weighted by the areas of each latitudinal strip, which vary as the cosine of latitude
Avg_insol = sum(pp.*cosd(lat))/sum(cosd(lat));

%Taking into account albedo: 
Avg_insol_absorbed = (1-albedo)*Avg_insol;
%Analytical calculation using the same So input as the insolation data: So = 1361 W/m^2
Avg_insol_absorbed_analytical = (1-albedo)*So/4;

yyaxis left
text(-75,175,['TOA insolation global annual average: ', sprintf('%5.2f',Avg_insol), ' W m^-^2'],'FontSize',12);
text(-75,160,['Absorbed TOA insolation = ', sprintf('%5.2f',Avg_insol_absorbed), ' W m^-^2'],'FontSize',12);
text(-75,145,['Analytical absorbed TOA insolation (S_o/4*(1-A)) = ', sprintf('%5.2f',Avg_insol_absorbed_analytical), ' W m^-^2'],'FontSize',12);