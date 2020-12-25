function out = dofyear2date_v2020(YYYY,first,varargin)
% Convert day of year to month/day and vice-versa
% Uses the datetime & duration classes
% Author: Dr. Tihomir Kostadinov, Nov. 7, 2020

%INPUTS: 
% Input either 1) YYYY - year (a scalar), and 2) a vector of days of year OR
% 1) YYYY - year (a scalar), 2) a vector of months, and 3) a corresponding vector of days

if ~isscalar(YYYY)
    error('Enter a single year in YYYY format');
end

if nargin==2
    %day of year given, convert to month, day
    dofyear = first(:);
    if ~isleap_vectorized(YYYY)
        qq = dofyear>365 | dofyear<1;
        if any(qq)
            warning('Invalid day of year detected');
            dofyear(qq) = NaN;
        end
    else
        qq = dofyear>366 | dofyear<1;
        if any(qq)
            warning('Invalid day of year detected');
            dofyear(qq) = NaN;
        end
    end
    dur = days(dofyear); %duration array
    t = datetime(YYYY,1,1,0,0,0)+dur-1;
    out = [t.Month, t.Day];
elseif nargin == 3
    %month and day given, convert to day of year
    mo = first(:);
    d = varargin{:}; d = d(:);
    
    c1 = (mo==1 | mo==3 | mo==5 | mo==7 | mo== 8 | mo ==10 | mo==12) & (d>=1 & d<=31);
    c2 = (mo==4 | mo==6 | mo==9 | mo==11 ) & (d>=1 & d<=30);
    
    if ~isleap_vectorized(YYYY)
        c3 = mo==2 & (d>=1 & d<=28);
    else
        c3 = mo==2 & (d>=1 & d<=29);
    end
    qq = c1 | c2 | c3;
   
    if any(~qq)
        warning('Invalid month/day combination detected');
        mo(~qq) = NaN;
        d(~qq) = NaN;
    end
    out = days(datetime(YYYY,mo,d,0,0,0) - datetime(YYYY,1,1,0,0,0) +1);
else
    error('Incorrect number of arguments');
end