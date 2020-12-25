function out = isleap_vectorized(years)
%Determines if a year is leap (output 1) or not (output 0) in the Gregorian calendar. 
%Input years - a vector of integer years. 
%Author: Dr. Tihomir Kostadinov

out =  ~(mod(years,4)~=0 | (mod(years, 100)==0 & mod(years,400)~=0));
