function [alpha_p_optim,A_p_optim,Annual_enegry_kWh] = solar_panel_orientation_optim(lat)
%Function to optimize the (fixed) orientation of a solar panel located at a given
%latitude. The annual power production is maximized using a combination of
%optimization techniques - namely, genetic algorithm, followed by fmincon
%using the 'interior-point' algorithm. 
%
%INPUT: 
%% lat - latitude of solar panel, degrees
%OUTPUTS:
%% alpha_p_optim - Optimal tilt of the solar panel, in degrees. Defined as angle between the normal to the plane of the solar panel and the local vertical direction.
%% A_p_optim - Optimal azimuth of the panel orientation, in degrees, clockwie from true North. This is the azimuth of the projection of the normal to the plane of the
%%       solar panel on the plane of the local horizon.
%% Annual_enegry_kWh - Annual energy produced at the solved optimal orientation. 
%%       This is the negative of the optimization cost function. Assumes a Pmax = 100 W solar panel at STC (1,000 W/m^2). 

%Author: Dr. Tihomir Kostadinov, Nov. 5-6, 2020

hybridopts = optimoptions('fmincon','OptimalityTolerance',1e-6);
options = optimoptions('ga','FunctionTolerance',1e-4,'Display','off',...
    'PopulationSize',200,'MaxStallGenerations', 20,'UseParallel',true, 'UseVectorized',false,...
    'HybridFcn',{'fmincon',hybridopts});
% Set 'UseParallel' to 'false' if you do not have the Parallel Computing Toolbox (c). 

if nargin==0 || isempty(lat)
    lat = 33+7/60+45/3600;
end

%Other defaults - choice of YYYY, TZ and lon, and Pmax should not influence the results much, if any.
lon = -117 -9/60 - 30/3600;
YYYY = 2021;
TZ = -8;
Pmax = 100;
dt = 5*60; %time step in seconds
to_plot = 0;

%Passing additional arguments via a nested function
f = @(YY)calc_cost_function_solar_production(YY);
[X, cost_fn_value] = ga(f,2,[],[],[],[],[0 0],[90 360],[],options);

alpha_p_optim = X(1);
A_p_optim = X(2);
Annual_enegry_kWh = -cost_fn_value;

%Nested cost function
    function cost_fn_value = calc_cost_function_solar_production(XX)
        %XX is a vector of the parameters to optimize for
        alpha_p = XX(1);
        A_p  = XX(2);       
        [~, Energy_annual_total] = solar_panel_annual(lat,lon,YYYY,TZ,alpha_p,A_p,Pmax,dt,to_plot);
        cost_fn_value = -Energy_annual_total;
    end %end of nested cost function
end