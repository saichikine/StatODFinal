% Saves station state data

%%
clear; close all; clc;

%% Loads

load('Project2a_Obs.mat');
load('Project2b_Obs.mat');

%%

tspan = Project2aObs(:,1); % 2a
%tspan = project2bObs(:,1); % 2b


rE = 6378.1363; % earth radius, km
spinRate = 7.29211585275553e-5; % earth spin rate, rad/s
theta0 = 0; %initial spin angle of Earth relative to ECI frame [degrees]
stationStatesLatLong = [rE+0.691750 -35.398333 148.981944; rE+0.834539 40.427222 355.749444; rE+1.07114904 35.247164 243.205;];

stationStates = zeros(6, 3, length(tspan));
stationOmega = [0, 0, spinRate]';
stationOmegas = repmat(stationOmega, 1, 1, length(stationStates));

%station inertial states
for i=1:3
   stationStates(1,i,:) = stationStatesLatLong(i,1)*sind(90 - stationStatesLatLong(i,2))*cos(spinRate*tspan + deg2rad(stationStatesLatLong(i,3) + theta0)); %station x
   stationStates(2,i,:) = stationStatesLatLong(i,1)*sind(90 - stationStatesLatLong(i,2))*sin(spinRate*tspan + deg2rad(stationStatesLatLong(i,3) + theta0)); %station y
   stationStates(3,i,:) = stationStatesLatLong(i,1)*cosd(90 - stationStatesLatLong(i,2)); %station z
   stationVel = cross(stationOmegas, stationStates(1:3,i,:), 1);
   stationStates(4:6,i,:) = stationVel;
end

% Save out station state data
save('stationStates.mat', 'stationStates');