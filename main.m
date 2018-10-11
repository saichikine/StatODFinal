% Sai Chikine
% ASEN 6080
% Homework 6 - main

%% 
clear; close all; clc;

%% Load Files

load('Project2_Prob2_truth_traj_50days.mat');
load('Project2a_Obs.mat');
load('measurementsPart2.mat');
load('stationStates2.mat');

%% Constants and stuff

AU =  149597870700/1000; %1 AU, [km]
c = 299792.458; %Speed of light, [km/s]

rE = 6378.1363; 

SRPFlux = 1357*1e6; %Solar pressure flux [W/km^2]
%SRPFlux = 0;

AoverM = 0.01/1e6; %SRP area to mass ratio [km^2/kg]

muSun = 132712440017.987; %[km^3/s^2]
muEarth = 3.986004415e5; %[km^3/s^2]

PhiSRP = SRPFlux*AU^2/1e6; 
PPhi = PhiSRP/c;

initialJD = 2456296.25; 

params = [muEarth, muSun, PPhi, AoverM, initialJD];

CRinit = 1;

% Truths for B-Plane Parameters
BPlaneParametersTruth = [14970.824, 9796.737];

%% Try integrating
tIntTest = [0, 50*86400];

opts = odeset('AbsTol', 1e-20, 'RelTol', 1e-13);

%ICs = [-274096790.0; -92859240.0; -40199490.0; 32.67; -8.94; -3.88; CRinit; reshape(eye(7), [], 1)];
ICs = [Xt_50(1,1:7)'; reshape(eye(7), [], 1)];

[t, X] = ode45(@(t,X) SRP3BDynamicsSTM(t, X, params), Tt_50, ICs, opts);

%% Plot satellite traj

[xSph, ySph, zSph] = sphere;
figure(1)
surf(rE*xSph,rE*ySph,rE*zSph)
hold on
plot3(X(:,1), X(:,2), X(:,3))
%scatter3(X0(1), X0(2), X0(3))
view([50 30])
xlabel('Inertial X, [km]')
ylabel('Inertial Y, [km]')
zlabel('Inertial Z, [km]')
title('Trajectory in ECI frame')

% %% Plot satellite traj
% 
% % [xSph, ySph, zSph] = sphere;
% % figure(2)
% % surf(rE*xSph,rE*ySph,rE*zSph)
% hold on
% plot3(Xt_50(:,1)-X(:,1), Xt_50(:,2)-X(:,2), Xt_50(:,3)-X(:,3))
% %scatter3(X0(1), X0(2), X0(3))
% view([50 30])
% xlabel('Inertial X, [km]')
% ylabel('Inertial Y, [km]')
% zlabel('Inertial Z, [km]')
% title('Trajectory in ECI frame DIFFERENCES')
% 
% %% Plot differences in position and vel
% figure
% subplot(211); hold on
% plot(t, Xt_50(:,1)-X(:,1), 'r-')
% plot(t, Xt_50(:,2)-X(:,2), 'g-')
% plot(t, Xt_50(:,3)-X(:,3), 'b-'); hold off
% 
% subplot(212); hold on
% plot(t, Xt_50(:,4)-X(:,4), 'r-')
% plot(t, Xt_50(:,5)-X(:,5), 'g-')
% plot(t, Xt_50(:,6)-X(:,6), 'b-'); hold off
% 
% %% Plot differences in diag elements of STM
% 
% figure
% hold on
% plot(t, Xt_50(:,8)-X(:,8))
% plot(t, Xt_50(:,16)-X(:,16))
% plot(t, Xt_50(:,24)-X(:,24))
% plot(t, Xt_50(:,32)-X(:,32))
% plot(t, Xt_50(:,40)-X(:,40))
% plot(t, Xt_50(:,48)-X(:,48))
% plot(t, Xt_50(:,56)-X(:,56))

%% UKF 

X02 = [-274096790.0; -92859240.0; -40199490.0; 32.67; -8.94; -3.88; 1.2];

P0 = diag([100 100 100 0.1 0.1 0.1 0.1]).^2;
R = diag([5/1000 0.5/1e6]).^2;
Q = 0*ones(7,7);

UKFInitStuff = {X02, P0, R, 60};
UKFNoiseStuff = {Q};
tic
UKF1 = UKF(UKFInitStuff, measurementsPart2(1:end,:), stationStates, UKFNoiseStuff, @SRP3BDynamics, @meas, [1, 2], params);
toc

%% Plot UKF Stuff

figure
subplot(211)
scatter(measurementsPart2(:,1), UKF1.postfits(1,:), 'DisplayName', 'Range Post-fits'); hold on
plot([0, measurementsPart2(end,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
plot([0, measurementsPart2(end,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
ylim([-0.03 0.03])
legend('Range Post-fits', '3\sigma Bounds')
ylabel('[km]')
title('Post-fits');
subplot(212)
scatter(measurementsPart2(:,1), UKF1.postfits(2,:)); hold on
plot([0, measurementsPart2(end,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--');
plot([0, measurementsPart2(end,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
ylim([-3e-6 3e-6])
legend('Range-Rate Post-fits', '3\sigma Bounds')
ylabel('[km/s]')
xlabel('Time [s]')

figure
subplot(211)
scatter(measurementsPart2(:,1), UKF1.prefits(1,:)); hold on
plot([0, measurementsPart2(end,1)], [3*5/1000, 3*5/1000], 'k--');
plot([0, measurementsPart2(end,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
ylim([-0.03 0.03])
legend('Range Post-fits', '3\sigma Bounds')
ylabel('[km]')
title('Pre-fits');
subplot(212)
scatter(measurementsPart2(:,1), UKF1.prefits(2,:)); hold on
plot([0, measurementsPart2(end,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--');
plot([0, measurementsPart2(end,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
ylim([-3e-6 3e-6])
legend('Range-Rate Post-fits', '3\sigma Bounds')
ylabel('[km/s]')
xlabel('Time [s]')
%% CKF

tic
CKFTest = CKFGeneral(X02, P0, zeros(7,1), R, measurementsPart2, stationStates, Q, 60, @SRP3BDynamicsSTM, @measPartials, @meas, params);
toc
%% Plot CKF Stuff

figure
hold on
plot(1:1:11064, CKFTest.Postfit_Residuals(1,:));
plot(1:1:11064, CKFTest.Postfit_Residuals(2,:)); hold off


%% Compare STM norms over 50 days
 
STM50DayNorms = zeros(1,length(t));

for i=1:length(t)
    STM50DayNorms(i) = norm(Xt_50(i,8:end)-X(i,8:end));
end

figure;
plot(t, STM50DayNorms)

%% Test iterated CKF
tic
CKFTestIterated = CKFPotterIterated(X02, P0, zeros(7,1), R, measurementsPart2(1:end,:), stationStates, @SRP3BDynamicsSTM, @measPartials, @meas, params, opts, 2, 0);
toc

%% 
figure
subplot(211)
scatter(measurementsPart2(:,1), CKFTestIterated.Prefit_Residuals(1,:,1)); hold on
plot([0, measurementsPart2(end,1)], [3*5/1000, 3*5/1000], 'k--');
plot([0, measurementsPart2(end,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
subplot(212)
scatter(measurementsPart2(:,1), CKFTestIterated.Prefit_Residuals(2,:,1)); hold on
plot([0, measurementsPart2(end,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--');
plot([0, measurementsPart2(end,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off

%% B-Plane Stuff (all data)

% % find B-plane intersection by dotting position vector with Shat vector
% 
% RSOI3 = 925000*3; % 3 times Earth sphere of influence distance, [km]
% RSOI3EventOpts = odeset('AbsTol', 1e-20, 'RelTol', 1e-13, 'Events', @(t,X) EventRSOI3(t,X,RSOI3));
% 
% ICs = [UKF1.stateEstimates(:,end); reshape(eye(7), [], 1)];
% 
% [tRSOI3, XRSOI3] = ode45(@(t,X) SRP3BDynamicsSTM(t, X, params), [measurementsPart2(end,1), inf], ICs, RSOI3EventOpts);

%% Plot distance

SOIIntegrationDistances = zeros(1,length(tRSOI3));
for i=1:length(tRSOI3)
    SOIIntegrationDistances(i) = norm(XRSOI3(i,1:3));
end

figure; hold on
scatter(tRSOI3, SOIIntegrationDistances, 'b.', 'DisplayName', 'Distance')
plot([tRSOI3(1) tRSOI3(end)], [RSOI3, RSOI3], 'DisplayName','3 RSOI Distance'); hold off
xlabel('Time past DCO [s]')
xlim([tRSOI3(1), tRSOI3(end)])
legend()

% %% Compute B Plane Things
% 
% BPlaneStuff = BPlaneCompute(XRSOI3(end,1:6)', muEarth);
% BPlaneDCM = BPlaneStuff{4};
% SHat = BPlaneStuff{1};
% 
% % Integrate to B-Plane crossing
% 
% BPlaneCrossingEventOpts = odeset('AbsTol', 1e-20, 'RelTol', 1e-13, 'Events', @(t,X) EventBPlaneCrossing(t,X,SHat));
% ICs = XRSOI3(end,:);
% 
% [tBPlane, XBPlane] = ode45(@(t,X) SRP3BDynamicsSTM(t, X, params), [tRSOI3(end), inf], ICs, BPlaneCrossingEventOpts);
% 
% %% Propagate final uncertainty to B
% P = UKF1.filterCovars(:,:,end);
% 
% STM = reshape(XBPlane(end,8:end),7,7);
% P = STM*P*STM';
% 
% %% Rotate Covars and pull out the on-plane parts
% PCrossingBPlaneFrame = BPlaneDCM*P(1:3,1:3)*BPlaneDCM';
% PCrossingBPlaneFrameTR = PCrossingBPlaneFrame(2:3,2:3);
%     
% %% B Plane Estimates
% BdotR = dot(BPlaneStuff{5}, BPlaneStuff{3});
% BdotT = dot(BPlaneStuff{5}, BPlaneStuff{2});
% 
% BPlaneParameters = [BdotR, BdotT];
% BPlaneErrors = BPlaneParameters - BPlaneParametersTruth;

%%

RSOI3 = 925000*3; % 3 times Earth sphere of influence distance, [km]

%% 50 Days Stuff

Index50Days = 3057;

tic
Results50Days = EstimateBPlaneParameters(measurementsPart2(Index50Days,1), UKF1.stateEstimates(:,Index50Days), UKF1.filterCovars(:,:,Index50Days), RSOI3, params);
toc
%% 100 Days Stuff

Index100Days = 5935;

tic
Results100Days = EstimateBPlaneParameters(measurementsPart2(Index100Days,1), UKF1.stateEstimates(:,Index100Days), UKF1.filterCovars(:,:,Index100Days), RSOI3, params);
toc

%% 150 Days Stuff

Index150Days = 8330;

tic
Results150Days = EstimateBPlaneParameters(measurementsPart2(Index150Days,1), UKF1.stateEstimates(:,Index150Days), UKF1.filterCovars(:,:,Index150Days), RSOI3, params);
toc
%% 200 Days Stuff

tic
Results200Days = EstimateBPlaneParameters(measurementsPart2(end,1), UKF1.stateEstimates(:,end), UKF1.filterCovars(:,:,end), RSOI3, params);
toc

%% Plot Ellipses

ellipse50Days = covarEllipse(Results50Days{2});
ellipse100Days = covarEllipse(Results100Days{2});
ellipse150Days = covarEllipse(Results150Days{2});
ellipse200Days = covarEllipse(Results200Days{2});

figure; hold on
plot(ellipse50Days(:,1)+Results50Days{1}(2), ellipse50Days(:,2)+Results50Days{1}(1), 'r-', 'DisplayName', '50 Days Ellipse');
plot(ellipse100Days(:,1)+Results100Days{1}(2), ellipse100Days(:,2)+Results100Days{1}(1), 'g-', 'DisplayName', '100 Days Ellipse')
plot(ellipse150Days(:,1)+Results150Days{1}(2), ellipse150Days(:,2)+Results150Days{1}(1), 'b-', 'DisplayName', '150 Days Ellipse')
plot(ellipse200Days(:,1)+Results200Days{1}(2), ellipse200Days(:,2)+Results200Days{1}(1), 'k-', 'DisplayName', '200 Days Ellipse')

scatter(Results50Days{1}(2), Results50Days{1}(1), 'rx', 'DisplayName', '50 Days Estimate')
scatter(Results100Days{1}(2), Results100Days{1}(1), 'gx', 'DisplayName', '100 Days Estimate')
scatter(Results150Days{1}(2), Results150Days{1}(1), 'bx', 'DisplayName', '150 Days Estimate')
scatter(Results200Days{1}(2), Results200Days{1}(1), 'kx', 'DisplayName', '200 Days Estimate')

xlabel('T')
ylabel('R')
legend()

%% Part 3 Test

%% Load
S = load('measurementsPart3.mat');
fieldNames = fieldnames(S);
fieldName = fieldNames{1};
measurementsPart3 = S.(fieldName);

load('stationStates3.mat');

%%

X03 = [-274096770.76544; -92859266.4499061; -40199493.6677441; 32.6704564599943; -8.93838913761049; -3.87881914050316; 1.0];
P0 = diag([100 100 100 0.1 0.1 0.1 0.1]).^2;
R = diag([5/1000 0.5/1e6]).^2;
Q = 0*ones(7,7);

UKFInitStuff_3 = {X03, P0, R};
UKFNoiseStuff = {Q};

filterIndex = 750;

tic
UKFP3 = UKF(UKFInitStuff_3, measurementsPart3(1:filterIndex,:), stationStates_3, UKFNoiseStuff, @SRP3BDynamics, @meas, [1, 2], params);
toc

%% Plot UKF Test for P3

% filterIndexStart = 13534;
% filterIndexEnd = 14200;

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.postfits(1,filterIndexStart:filterIndexEnd,1), '.','DisplayName', 'Range Post-fits, Station 1');
r2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.postfits(1,filterIndexStart:filterIndexEnd,2), '.','DisplayName', 'Range Post-fits, Station 2');
r3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.postfits(1,filterIndexStart:filterIndexEnd,3), '.','DisplayName', 'Range Post-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-100 150])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Post-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.postfits(2,filterIndexStart:filterIndexEnd,1), '.','DisplayName', 'Range-Rate Post-fits, Station 1');
rr2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.postfits(2,filterIndexStart:filterIndexEnd,2), '.','DisplayName', 'Range-Rate Post-fits, Station 2');
rr3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.postfits(2,filterIndexStart:filterIndexEnd,3), '.','DisplayName', 'Range-Rate Post-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-10e-4 5e-4])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.prefits(1,filterIndexStart:filterIndexEnd,1), '.','DisplayName', 'Range Pre-fits, Station 1');
r2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.prefits(1,filterIndexStart:filterIndexEnd,2), '.','DisplayName', 'Range Pre-fits, Station 2');
r3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.prefits(1,filterIndexStart:filterIndexEnd,3), '.','DisplayName', 'Range Pre-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-100 150])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Pre-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.prefits(2,filterIndexStart:filterIndexEnd,1), '.','DisplayName', 'Range-Rate Pre-fits, Station 1');
rr2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.prefits(2,filterIndexStart:filterIndexEnd,2), '.','DisplayName', 'Range-Rate Pre-fits, Station 2');
rr3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3.prefits(2,filterIndexStart:filterIndexEnd,3), '.','DisplayName', 'Range-Rate Pre-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-10e-4 5e-4])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

%% Try UKF with DMC

tau = 60; %time constant for DMC accelerations [s]
paramsDMC = [muEarth, muSun, PPhi, AoverM, initialJD, tau];

filterIndexStart = 13534;
filterIndexEnd = 14036;

X03 = [UKFP3.stateEstimates(1:6,filterIndexStart); zeros(3,1); UKFP3.stateEstimates(7,filterIndexStart)];
P0 = diag([100 100 100 0.1 0.1 0.1 1e-12 1e-12 1e-12 0.1]).^2;
R = diag([5/1000 0.5/1e6]).^2;
Q = diag([0 0 0]);
UKFInitStuff_3 = {X03, P0, R, filterIndexStart, filterIndexEnd};
UKFNoiseStuff = {Q};

tic
UKFP3TestDMC = UKF_DMC(UKFInitStuff_3, measurementsPart3(filterIndexStart:filterIndexEnd,:), stationStates_3(:,:,filterIndexStart:filterIndexEnd), UKFNoiseStuff, @SRP3BDynamicsDMC, @meas, [1, 2], paramsDMC);
toc

%% UKF w DMC Test Plots

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.postfits(1,:,1), '.','DisplayName', 'Range Post-fits, Station 1');
r2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.postfits(1,:,2), '.','DisplayName', 'Range Post-fits, Station 2');
r3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.postfits(1,:,3), '.','DisplayName', 'Range Post-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-10 50])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Post-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.postfits(2,:,1), '.','DisplayName', 'Range-Rate Post-fits, Station 1');
rr2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.postfits(2,:,2), '.','DisplayName', 'Range-Rate Post-fits, Station 2');
rr3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.postfits(2,:,3), '.','DisplayName', 'Range-Rate Post-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-3e-3 3e-3])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.prefits(1,:,1), '.','DisplayName', 'Range Pre-fits, Station 1');
r2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.prefits(1,:,2), '.','DisplayName', 'Range Pre-fits, Station 2');
r3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.prefits(1,:,3), '.','DisplayName', 'Range Pre-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-10 50])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Pre-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.prefits(2,:,1), '.','DisplayName', 'Range-Rate Pre-fits, Station 1');
rr2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.prefits(2,:,2), '.','DisplayName', 'Range-Rate Pre-fits, Station 2');
rr3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFP3TestDMC.prefits(2,:,3), '.','DisplayName', 'Range-Rate Pre-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-3e-3 3e-3])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

%% UKF after DMC

X0Post = [UKFP3TestDMC.stateEstimates(1:6,end); UKFP3TestDMC.stateEstimates(10,end)];
P0 = diag([100 100 100 0.1 0.1 0.1 0.1]).^2;
R = diag([5/1000 0.5/1e6]).^2;
Q = 0*ones(7,7);

filterIndexStart = 14036;
filterIndexEnd = length(measurementsPart3);

UKFInitStuff_3 = {X0Post, P0, R, filterIndexStart, filterIndexEnd};
UKFNoiseStuff = {Q};

tic
UKFPost = UKF(UKFInitStuff_3, measurementsPart3(filterIndexStart:filterIndexEnd,:), stationStates_3(:,:,filterIndexStart:filterIndexEnd), UKFNoiseStuff, @SRP3BDynamics, @meas, [1, 2], params);
toc

%% UKF Post DMC Plots

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.postfits(1,:,1), '.','DisplayName', 'Range Post-fits, Station 1');
r2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.postfits(1,:,2), '.','DisplayName', 'Range Post-fits, Station 2');
r3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.postfits(1,:,3), '.','DisplayName', 'Range Post-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-10 50])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Post-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.postfits(2,:,1), '.','DisplayName', 'Range-Rate Post-fits, Station 1');
rr2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.postfits(2,:,2), '.','DisplayName', 'Range-Rate Post-fits, Station 2');
rr3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.postfits(2,:,3), '.','DisplayName', 'Range-Rate Post-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-3e-3 3e-3])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.prefits(1,:,1), '.','DisplayName', 'Range Pre-fits, Station 1');
r2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.prefits(1,:,2), '.','DisplayName', 'Range Pre-fits, Station 2');
r3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.prefits(1,:,3), '.','DisplayName', 'Range Pre-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-10 50])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Pre-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.prefits(2,:,1), '.','DisplayName', 'Range-Rate Pre-fits, Station 1');
rr2 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.prefits(2,:,2), '.','DisplayName', 'Range-Rate Pre-fits, Station 2');
rr3 = scatter(measurementsPart3(filterIndexStart:filterIndexEnd,1), UKFPost.prefits(2,:,3), '.','DisplayName', 'Range-Rate Pre-fits, Station 3');
b1 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(filterIndexStart,1), measurementsPart3(filterIndexEnd,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-3e-3 3e-3])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

%% Combined

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(:,1), [UKFP3.postfits(1,1:13533,1), UKFP3TestDMC.postfits(1,2:end,1), UKFPost.postfits(1,:,1)], '.','DisplayName', 'Range Post-fits, Station 1');
r2 = scatter(measurementsPart3(:,1), [UKFP3.postfits(1,1:13533,2), UKFP3TestDMC.postfits(1,2:end,2), UKFPost.postfits(1,:,2)], '.','DisplayName', 'Range Post-fits, Station 2');
r3 = scatter(measurementsPart3(:,1), [UKFP3.postfits(1,1:13533,3), UKFP3TestDMC.postfits(1,2:end,3), UKFPost.postfits(1,:,3)], '.','DisplayName', 'Range Post-fits, Station 1');
b1 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-10 50])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Post-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(:,1), [UKFP3.postfits(2,1:13533,1), UKFP3TestDMC.postfits(2,2:end,1), UKFPost.postfits(2,:,1)], '.','DisplayName', 'Range-Rate Post-fits, Station 1');
rr2 = scatter(measurementsPart3(:,1), [UKFP3.postfits(2,1:13533,2), UKFP3TestDMC.postfits(2,2:end,2), UKFPost.postfits(2,:,2)], '.','DisplayName', 'Range-Rate Post-fits, Station 2');
rr3 = scatter(measurementsPart3(:,1), [UKFP3.postfits(2,1:13533,3), UKFP3TestDMC.postfits(2,2:end,3), UKFPost.postfits(2,:,3)], '.','DisplayName', 'Range-Rate Post-fits, Station 1');
b1 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-3e-3 3e-3])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

figure
subplot(211); hold on
r1 = scatter(measurementsPart3(:,1), [UKFP3.prefits(1,1:13533,1), UKFP3TestDMC.prefits(1,2:end,1), UKFPost.prefits(1,:,1)], '.','DisplayName', 'Range Pre-fits, Station 1');
r2 = scatter(measurementsPart3(:,1), [UKFP3.prefits(1,1:13533,2), UKFP3TestDMC.prefits(1,2:end,2), UKFPost.prefits(1,:,2)], '.','DisplayName', 'Range Pre-fits, Station 2');
r3 = scatter(measurementsPart3(:,1), [UKFP3.prefits(1,1:13533,3), UKFP3TestDMC.prefits(1,2:end,3), UKFPost.prefits(1,:,3)], '.','DisplayName', 'Range Pre-fits, Station 1');
b1 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [3*5/1000, 3*5/1000], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [-3*5/1000, -3*5/1000], 'k--'); hold off
%ylim([-10 50])
legend([r1 r2 r3 b1])
ylabel('[km]')
title('Pre-fits');
subplot(212); hold on
rr1 = scatter(measurementsPart3(:,1), [UKFP3.prefits(2,1:13533,1), UKFP3TestDMC.prefits(2,2:end,1), UKFPost.prefits(2,:,1)], '.','DisplayName', 'Range-Rate Pre-fits, Station 1');
rr2 = scatter(measurementsPart3(:,1), [UKFP3.prefits(2,1:13533,2), UKFP3TestDMC.prefits(2,2:end,2), UKFPost.prefits(2,:,2)], '.','DisplayName', 'Range-Rate Pre-fits, Station 2');
rr3 = scatter(measurementsPart3(:,1), [UKFP3.prefits(2,1:13533,3), UKFP3TestDMC.prefits(2,2:end,3), UKFPost.prefits(2,:,3)], '.','DisplayName', 'Range-Rate Pre-fits, Station 1');
b1 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [3*0.5/1e6, 3*0.5/1e6], 'k--', 'DisplayName', '3\sigma Bounds');
b2 = plot([measurementsPart3(1,1), measurementsPart3(end,1)], [-3*0.5/1e6, -3*0.5/1e6], 'k--'); hold off
%ylim([-3e-3 3e-3])
legend([rr1 rr2 rr3 b1])
ylabel('[km/s]')
xlabel('Time [s]')

%% Try estimating station positions

X03s = [-274096770.76544; -92859266.4499061; -40199493.6677441; 32.6704564599943; -8.93838913761049; -3.87881914050316; ...
    stationStates_3(1:3,1,1); stationStates_3(1:3,2,1); stationStates_3(1:3,3,1); 1.0];
P0s = diag([100 100 100 0.1 0.1 0.1 1e-12*ones(1,3) 10*ones(1,3) 1e-8*ones(1,3) 0.1]).^2;
R = diag([5/1000 0.5/1e6]).^2;
Q = 0*ones(length(X03s),length(X03s));

UKFInitStuff_3 = {X03s, P0s, R};
UKFNoiseStuff = {Q};

filterIndexStart = 1;
filterIndexEnd = 750;

tic
UKF3s = UKFs(UKFInitStuff_3, measurementsPart3(filterIndexStart:filterIndexEnd,:), UKFNoiseStuff, @SRP3BDynamics_s, @meas_s, [1, 2], params);
toc

%% Station Estimation Results

plotFilterResiduals(measurementsPart3, UKF3s.prefits, UKF3s.postfits, [filterIndexStart filterIndexEnd])

