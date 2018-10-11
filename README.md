# StatODFinal
Final project for ASEN 6020 Statistical Orbit Determination. Uses Unscented Kalman Filter with Dynamic Model Compensation to estimate spacecraft state with unknown perturbations. Known dynamics include circular restricted three-body problem equations of motion and cannonball model solar radiation pressure. Also uses a classical Kalman filter (CKF) on more well behaved data.

See file `UKF_DMCParams.m` for unscented Kalman filter code. See `CKFGeneral.m` for classical Kalman filter. 
