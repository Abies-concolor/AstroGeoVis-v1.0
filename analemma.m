%Plot the analemma and related quantities
%Author: Tihomir Kostadinov, 2014 - 2020 (major re-write and improvements in 2020)
%
%%Inputs
YYYY = 2021; %Year chosen for the calculations. 
lon = -117.25; %CSU San Marcos longitude in degrees (Negative in Western Hemisphere)
TZ = -8; %Time zone offset from Universal Time in hours (positibve East of Greenwich; add 1 hour for Daylight Saving Time)

%Constant
ks = 1.00273790935; %This is the scaling between an interval of mean solar(civil) time and an interval of sidereal time. 
%In other words, ds/dt = ks where s is the sidereal time
%Note that this value differs slighlty in Meeus's sidereal time algortihm data, and it varies slowly with time, 
%for details see Explanatory Supplement to the Astronomical Almanac, 3rd Ed., pg. 81

t = datetime(YYYY,1,1,0,0,0):1/86400:datetime(YYYY,12,31,23,59,59);
t_local = t +TZ/24;
JD = date2jd_vec(t.Year,t.Month,t.Day,t.Hour,t.Minute,t.Second,'G'); %Gregorian calendar assumed

[RA_Sun, dec_Sun, Rvector] = solar_coord(JD); %RA & dec in degrees, Rvector in AU
[local_s, So] = sidereal_time(JD,lon); %sidereal time in hours

%Hour angle of the true Sun:
HA_true_Sun = local_s - RA_Sun/15; 
HA_true_Sun = mod(HA_true_Sun,24); %bring to a range of 0 to 24 hrs. 

%%%%%%% Explore the rate of change of solar right ascension (RA) %%%%%%%%%%%%%
alpha = RA_Sun;
dalpha_dt = mod(diff(alpha),360); %Takes care of circular discontinuity; units are degrees/sec as dt = 1 s in t sampling
dalpha_dt_min_per_day = 86400*(dalpha_dt/15)*60; %in minutes (of time, not arcminutes of angle) of RA per day

%figure
%plot(t(2:end),dalpha_dt_min_per_day)

%(1)
% H = s - alpha
% Thus dH/dt = ds/dt - dalpha/dt, and ds/dt = 1.002737... = ks above
% (2)
%By definition, EqT = HA_true_Sun - HA_Mean_Sun
%Differentiate the above and get: dEqT/dt = HA_true/dt - dHA_mean/dt, and
%dHA_mean/dt = 1 by definition, so we substitude the differential of (1) above & get: 
%dEqT/dt = (ks-1) - dRA_Sun/dt -- see plot and interpretation below:

dEqT_dt = (ks-1)*86400/60 - dalpha_dt_min_per_day; %ks and '1' have units of second/second, which are converted to min/day
dsdt_min_day = (ks-1)*86400/60; 
%The rate of change of the equation of time would be equal to exactly 0 IF 
%dRA_Sun/dt is equal to (ks-1),i.e. if the Sun's right ascension changed
%uniformly. It does not due to the reasons for the analemma - eccentricity and obliquity!

figure
yyaxis left
plot(t(2:end),dalpha_dt_min_per_day,'LineWidth',1.5)
hold on
plot([t(2) t(end)],[dsdt_min_day dsdt_min_day],'-.')
ylim([dsdt_min_day-0.5 dsdt_min_day+0.5])
ylabel('d\alpha\_true\_Sun/dt, min/day')
yyaxis right
plot(t(2:end),dEqT_dt,'LineWidth',1.5)
ylim([-0.5 0.5])
yticks(-0.5:.1:0.5)
ylabel('dEqT/dt, min/day')
set(gca,'XTick',datetime(YYYY,1:12,15,12,0,0))
xtickformat('MMM.dd')
xtickangle(45)
xlabel('Date')
legend({'d\alpha\_true\_Sun/dt, min/day','Mean Solar day - Sidereal day, min.','dEqT/dt, min/day'})
title(['Time Derivatives of the true Sun \alpha and the Equation of Time'])

%Compute and plot the Equation of Time and the Analemma
longitude_correction = TZ - lon/15; %offset from central longitude of time zone, i.e. the constant offset between the 
%time zone's central meridian (typically multiple of 15 degrees = 1 hour) and the actual observer's longitude.
HA_mean_Sun = mod((t_local.Hour + t_local.Minute/60 + t_local.Second/3600 - 12 -longitude_correction)',24); 
%Relationship of local civilian time to hour angle of the mean Sun
%Local Hour angle of the Mean Sun + 12 + longitude_correction = local
%civilian time; t_local represents the local civilian time, the above is just this
%relationship solved for HA_mean_Sum

EqT = HA_true_Sun - HA_mean_Sun; %By definition, this is the equation of time. 
EqT(EqT>1) = EqT(EqT>1)-24;
EqT(EqT<-1) = EqT(EqT<-1)+24;
%Therefore one way to define The Equation of Time is: It is equal to the hour angle of the real Sun at local MEAN solar
%noon (when the hour angle of the MEAN Sun is 0 hours). 
%Students can be tasked to reproduce the Equation of time using this definiton.
%Local mean solar noon occurs at t = 12 h + longitude_correction

%Plot the Equation of Time
figure
plot(t,EqT*60,'LineWidth',1.5)
set(gca,'XTick',datetime(YYYY,1:12,15,12,0,0))
hold on
plot([t(1) t(end)],[0 0],'-.','LineWidth',1.5)
ylabel('Equation of Time, minutes')
xtickformat('MMM.dd')
xlabel('Date')
xtickangle(45)
yticks(-15:5:15)
title(['The Equation of Time for ', num2str(YYYY)])

%Plot the analemma
figure
plot(EqT*60,dec_Sun,'LineWidth',1.5)
hold on
grid on
xlim([-15 17])
xticks([-15:5:15])
ylim([-25 25])
yticks([-25:5:25])
plot([-15 17],[0 0],'k-.','LineWidth',1.5);
plot([0 0],[-25 25],'k-.','LineWidth',1.5);
grid on
xlabel('Equation of Time = True Sun HA - Mean Sun HA, min.')
ylabel('Solar declination, deg')
%Plot the beginning of each month in red 
mobegin_idx = find(ismember(t,datetime(YYYY,1:12,1,0,0,0)))'; 
hold on
plot(EqT(mobegin_idx)*60,dec_Sun(mobegin_idx),'ro')
textx = EqT(mobegin_idx)*60;
texty = dec_Sun(mobegin_idx);
textx(11) = textx(11) - 5;
textx(1) = textx(1) - 5;
textx(7) = textx(7) - 5;
str = {'Jan 1 \rightarrow',' \leftarrow Feb 1',' \leftarrow Mar 1',' \leftarrow Apr 1',...
    ' \leftarrow May 1',' \leftarrow Jun 1','Jul 1 \rightarrow',' \leftarrow Aug 1',...
    ' \leftarrow Sep 1',' \leftarrow Oct 1','Nov 1 \rightarrow',' \leftarrow Dec 1'};
text(textx,texty,str,'FontSize',14);
title(['Analemma for the Year ',num2str(YYYY)]);
set(gca,'FontSize',14)

%Plot the analemma on the same x and y axis scale - both axes in degrees
figure
plot(EqT*15,dec_Sun)
set(gca,'DataAspectRatio', [ 1 1 1 ]);
xlim([-5 5])
xticks([-5:2.5:5])
ylim([-25 25])
yticks([-25:5:25])
hold on
plot([-5 5],[0 0],'k-.');
plot([0 0],[-25 25],'k-.');
grid on
xlabel('Equation of Time = True Sun HA - Mean Sun HA, deg.')
ylabel('Solar declination, deg')
xticklabels(sprintf('%3.1f\n',[-5:2.5:5]))
xtickangle(90)
title(['High Resolution Analemma for the year ',num2str(YYYY)]);