function JD = date2jd_vec(Y,M,D,h,m,s,type)
%Based on Meeus, J. (1998), Astronomical Algorithms, Ch. 7
%
% INPUT - The numeric inputs be Nx1 column vectors, where N is the number of samples to be converted
% Y:year
% M:month - 1 - Jan, 2 - Feb, etc
% D: day of the month
% h,m,s : hour, minute, seconds, in 24 hr system;
% The above are assumed to be given in Universal Time, apply any time zone offset before passing to this function!!!
% type: 1x1 character array or, 'J' --> Julian; 'G'-->Gregorian calendar
%
% OUTPUT:
% JD: The Julian Day number, days since the beginning of -4712, January 1, Noon UT
% -4712 corresponds to 4713 BC, see Meeus Ch. 7 - astronomical vs historical counting
%
% Caution: Input checking is basic, not guaramteed to be rigorous and/or complete.
% Written in MATLAB(r), vectorization & input checking by Dr. Tihomir S. Kostadinov

%Basic input checking
Y = Y(:);
M = M(:);
D = D(:);
h = h(:);
m = m(:);
s = s(:);

if ~isscalar(unique([numel(Y) numel(M) numel(D) numel(h) numel(m) numel(s)]))
    error('Invalid input')
end

%May be unnecessarily strict, but good practice
if any(Y<-4712) | any(M<1 | M>12) | any(D<0 | D>31) | any(h<0 | h>24) | any(m<0 | m> 60) |any (s<0 | s>60)
    error('Double check that all date components are input correctly and have valid ranges.')
end
%End input checking

D  = D + (h + m/60 + s/3600)/24; %convert to decimal day
%Day can also be supplied as a fraction, then h,m,s should be set to 0's;

qq = M==1 | M==2;
Y(qq) = Y(qq)-1;
M(qq) = M(qq)+12;

if strncmpi(type,'G',1)
    A = floor(Y/100);
    B = 2-A+floor(A/4);
elseif strncmpi(type,'J',1)
    B = 0;
else
    error('Calendar type not recognized, specify the type of Calendar, Gregorian(''G'') or Julian(''J'')');
end

JD = floor(365.25*(Y+4716)) + floor(30.6001*(M+1)) + D + B - 1524.5;