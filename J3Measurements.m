% Measurement Simulation

%%
clear all; close all; clc;

%% Constants
%Keplerian orbital elements
a = 10000; %semi-major axis, [km]
e = 0.001; %eccentricity
i = deg2rad(40); %inclination angle, [radians]
RAAN = deg2rad(80); %right ascension of ascending node, [radians]
w = deg2rad(40); %argument of periapsis [radians]
nu = deg2rad(0); %initial true anomaly [radians]

mu = 3.986004415e5; %Earth gravitational parameter, [km^4/s^2]
rE = 6378.136; %Earth radius, [km]
J2 = 0.0010826269; %Earth J2 parameter
J3 = -0.0000025323; %Earth j3 parameter
numOrbits = 15; %number of orbits to integrate
orbitPeriod = 2*pi*sqrt(a^3/mu); %orbital period, [s]
simTime = numOrbits*orbitPeriod; %total simulation time, [s]
tspan = 0:10:simTime; %10 second steps for measurements

spinRate = 7.2921158553e-5; %Earth spin rate, [rad/s]

params = [mu, J2, rE, J3];

%% Convert initial Keplerian constants to initial Cartesian state
p = a*(1-e^2); %semi-latus rectum, [km]
r = p/(1+e*cos(nu)); %orbit radius, [km]
rp = a*(1-e); %periapsis radius, [km]

h = sqrt(mu*a*(1-e^2)); %angular momentum

x = r*(cos(RAAN)*cos(w+nu) - sin(RAAN)*sin(w+nu)*cos(i)); %x-position, [km]
y = r*(sin(RAAN)*cos(w+nu) + cos(RAAN)*sin(w+nu)*cos(i)); %y-position, [km]
z = r*(sin(i)*sin(w+nu)); %z-position, [km]

xdot = x*h*e/(r*p)*sin(nu) - h/r*(cos(RAAN)*sin(w+nu) + sin(RAAN)*cos(w+nu)*cos(i)); %x-velocity, [km/s]
ydot = y*h*e/(r*p)*sin(nu) - h/r*(sin(RAAN)*sin(w+nu) - cos(RAAN)*cos(w+nu)*cos(i)); %y-velocity, [km/s]
zdot = z*h*e/(r*p)*sin(nu) + h/r*sin(i)*cos(w+nu);

X0 = [x y z xdot ydot zdot]'; %initial state vector
%% Compute reference trajectory (no perturbation)

%deltx0 = [1 0 0 0 10/1000 0]'; %[r v]'
deltx0 = [0 0 0 0 0 0 ]';
X0pert = X0 + deltx0;

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[t, XJ3] = ode45(@(t,X) J3Dynamics(t, X, params), tspan, X0pert, options);


%% Plot trajectory

[xSph, ySph, zSph] = sphere;
figure(1)
surf(rE*xSph,rE*ySph,rE*zSph)
hold on
plot3(XJ3(:,1), XJ3(:,2), XJ3(:,3))
scatter3(X0(1), X0(2), X0(3))
view([50 30])
xlabel('Inertial X, [km]')
ylabel('Inertial Y, [km]')
zlabel('Inertial Z, [km]')
title('Trajectory in ECI frame')

%% Generate station states

theta0 = 122; %initial spin angle of Earth relative to ECI frame [degrees]
stationStatesLatLong = [rE -35.39833 148.981944; rE 40.427222 355.74944; rE 35.247164 243.205;];

% stationECEFPositions = [
%  rE*cosd(stationStatesLatLong(1,2)*cosd(stationStatesLatLong(1,3) + theta), rE*cosd(stationStatesLatLong(1,2)*sind(stationStatesLatLong(1,3) + theta), rE*sind(stationStatesLatLong(1,2);
%  rE*cosd(stationStatesLatLong(2,2)*cosd(stationStatesLatLong(2,3) + theta), rE*cosd(stationStatesLatLong(2,2)*sind(stationStatesLatLong(2,3) + theta), rE*sind(stationStatesLatLong(2,2);
%  rE*cosd(stationStatesLatLong(3,2)*cosd(stationStatesLatLong(3,3) + theta), rE*cosd(stationStatesLatLong(3,2)*sind(stationStatesLatLong(3,3) + theta), rE*sind(stationStatesLatLong(3,2)];

stationStates = zeros(6, 3, length(tspan));
stationOmega = [0, 0, spinRate]';
stationOmegas = repmat(stationOmega, 1, 1, length(stationStates));

%station inertial states
for i=1:3
   stationStates(1,i,:) = rE*sind(90 - stationStatesLatLong(i,2))*cos(spinRate*tspan + deg2rad(stationStatesLatLong(i,3) + theta0)); %station x
   stationStates(2,i,:) = rE*sind(90 - stationStatesLatLong(i,2))*sin(spinRate*tspan + deg2rad(stationStatesLatLong(i,3) + theta0)); %station y
   stationStates(3,i,:) = rE*cosd(90 - stationStatesLatLong(i,2)); %station z
   stationVel = cross(stationOmegas, stationStates(1:3,i,:), 1);
   stationStates(4:6,i,:) = stationVel;
%    stationStates(4,i,:) = -spinRate*rE*sind(90 - stationStatesLatLong(i,2))*sin(spinRate*tspan + deg2rad(stationStatesLatLong(i,3) + theta0)); %station xdot
%    stationStates(5,i,:) = spinRate*rE*sind(90 - stationStatesLatLong(i,2))*cos(spinRate*tspan + deg2rad(stationStatesLatLong(i,3) + theta0)); %station xdot
%    stationStates(6,i,:) = 0; %station zdot
end

%% Simulate measurements

elevationMask = 10; %station elevation mask, stations can only see spacecraft if it is at least 10 degrees above station horizon [degrees]

y = zeros(3, 2, length(tspan));
stationElevationAngles = zeros(3, 1, length(tspan));

meas_Aformat = [];

% Noise characteristics for measurements
sigmaRange = 0.001; % 1m [km]
sigmaRangeRate = 1e-6; % 1mm/s [km/s]

for i = 1:3
    for j = 1:length(tspan)
        %RStationSC projected onto RStation >= RStationSC*sin(10 degrees) for
        %visibility <-- doesn't work, not sure why
        
        rangeNoise = normrnd(0, sigmaRange);
        rangeRateNoise = normrnd(0, sigmaRangeRate);
        
        RStationSC = XJ3(j,1:3)' - stationStates(1:3,i,j); %vector from station to space station
        RStationUnit = stationStates(1:3,i,j)/norm(stationStates(1:3,i,j)); %space station position unit vector
        
        stationElevationAngles(i,1,j) = 90 - atan2d(norm(cross(RStationSC, stationStates(1:3,i,j))), dot(RStationSC, stationStates(1:3,i,j)));
        
        if stationElevationAngles(i,1,j) >= 10
            y(i,1,j) = norm(RStationSC) + rangeNoise; %range measurement
            y(i,2,j) = dot(RStationSC',(XJ3(j,4:6)' - stationStates(4:6,i,j)))/norm(RStationSC) + rangeRateNoise; %range rate measurement
            
            % Save Andrew format data
            meas_Aformat = [meas_Aformat;
                            tspan(j), i, norm(RStationSC) + rangeNoise, ...
                            ((RStationSC'*(XJ3(j,4:6)' - stationStates(4:6,i,j)))/norm(RStationSC)) + rangeRateNoise];
        else
            y(i,:,j) = NaN; %no measurement
        end
    end
end

%% Save out data

%Save out measurements with noise for hw3

sigmarange = 0.001; %sigma range is 1m
sigmarangerate = 1/(1e6); %sigma range rate is 1 mm/s

youtJ3 = zeros(size(y));
youtJ3(:,1,:) = y(:,1,:) + normrnd(0, sigmarange, 3, 1, length(tspan));
youtJ3(:,2,:) = y(:,2,:) + normrnd(0, sigmarangerate, 3, 1, length(tspan));

save('measurementsJ3_saiformat.mat', 'youtJ3');

% Convert measurements to Andrew's format
measAFormatSorted = sortrows(meas_Aformat,1);
save('measurements_A.mat', 'measAFormatSorted');

% Save out station state data
save('stationStates.mat', 'stationStates');

%Save out true orbit data for hw3
save('trueorbitdataJ3.mat', 'XJ3');

%% Plot measurements

figure
subplot(2,3,1)
plot(tspan, reshape(y(1,1,:), 1, []))
xlabel('Time [s]')
ylabel('Range [km]')
title('Station 1 Range')

subplot(2,3,2)
plot(tspan, reshape(y(2,1,:), 1, []))
xlabel('Time [s]')
ylabel('Range [km]')
title('Station 2 Range')

subplot(2,3,3)
plot(tspan, reshape(y(3,1,:), 1, []))
xlabel('Time [s]')
ylabel('Range [km]')
title('Station 3 Range')

subplot(2,3,4)
plot(tspan, reshape(y(1,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 1 Range Rate')

subplot(2,3,5)
plot(tspan, reshape(y(2,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 2 Range Rate')

subplot(2,3,6)
plot(tspan, reshape(y(3,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 3 Range Rate')

% Plot elevation angles
figure
subplot(3,1,1)
hold on
plot(tspan, reshape(stationElevationAngles(1,1,:), 1, []))
plot([0, tspan(end)], [10, 10])
xlabel('Time [s]')
ylabel('Elevation Angle [degrees]')
title('Station 1 Elevation')
hold off

subplot(3,1,2)
hold on
plot(tspan, reshape(stationElevationAngles(2,1,:), 1, []))
plot([0, tspan(end)], [10, 10])
xlabel('Time [s]')
ylabel('Elevation Angle [degrees]')
title('Station 2 Elevation')
hold off

subplot(3,1,3)
hold on
plot(tspan, reshape(stationElevationAngles(3,1,:), 1, []))
plot([0, tspan(end)], [10, 10])
xlabel('Time [s]')
ylabel('Elevation Angle [degrees]')
title('Station 3 Elevation')
hold off

%% Real units

fref = 8.44; %reference transmit frequency [GHz]
c = 299792458/1000; %speed of light [km/s]
 
freqData = -2*y(:,2,:)/c*fref; %range rate converted to doppler shift freq data

rangeUnits = 221/749*y(:,1,:)/c*fref; %range converted to range units

%Plot new data

figure
subplot(2,3,1)
plot(tspan, reshape(rangeUnits(1,1,:), 1, []))
xlabel('Time [s]')
ylabel('Range [RU]')
title('Station 1 Range')

subplot(2,3,2)
plot(tspan, reshape(rangeUnits(2,1,:), 1, []))
xlabel('Time [s]')
ylabel('Range [RU]')
title('Station 2 Range')

subplot(2,3,3)
plot(tspan, reshape(rangeUnits(3,1,:), 1, []))
xlabel('Time [s]')
ylabel('Range [RU]')
title('Station 3 Range')

subplot(2,3,4)
plot(tspan, reshape(freqData(1,1,:), 1, []))
xlabel('Time [s]')
ylabel('Doppler Shift [GHz]')
title('Station 1 Range Rate')

subplot(2,3,5)
plot(tspan, reshape(freqData(2,1,:), 1, []))
xlabel('Time [s]')
ylabel('Doppler Shift [GHz]')
title('Station 2 Range Rate')

subplot(2,3,6)
plot(tspan, reshape(freqData(3,1,:), 1, []))
xlabel('Time [s]')
ylabel('Doppler Shift [GHz]')
title('Station 3 Range Rate')

%% Data with noise

sigma = 0.5/(1e6);

rangeRateNoise = normrnd(0, sigma, 3, 1, length(tspan));
noisyRangeRate = y(:,2,:) + rangeRateNoise;

% Plot noisy range rate data

figure
subplot(3,1,1)
hold on
plot(tspan, reshape(noisyRangeRate(1,1,:), 1, []))
plot(tspan, reshape(y(1,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 1 Range Rate')
legend('Noisy Range Rate Data', 'Original Range Rate Data')
hold off

subplot(3,1,2)
hold on
plot(tspan, reshape(noisyRangeRate(2,1,:), 1, []))
plot(tspan, reshape(y(2,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 2 Range Rate')
hold off

subplot(3,1,3)
hold on
plot(tspan, reshape(noisyRangeRate(3,1,:), 1, []))
plot(tspan, reshape(y(3,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 3 Range Rate')
hold off

%% Noisy minus original measurements
figure
subplot(3,1,1)
hold on
plot(tspan, reshape(noisyRangeRate(1,1,:)-y(1,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 1 Range Rate Difference')
legend('Noisy minus original range rates')
hold off

subplot(3,1,2)
hold on
plot(tspan, reshape(noisyRangeRate(2,1,:)-y(2,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 2 Range Rate Difference')
hold off

subplot(3,1,3)
hold on
plot(tspan, reshape(noisyRangeRate(3,1,:)-y(3,2,:), 1, []))
xlabel('Time [s]')
ylabel('Range Rate [km/s]')
title('Station 3 Range Rate Difference')
hold off