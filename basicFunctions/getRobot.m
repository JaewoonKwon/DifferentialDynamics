function initialRobot = getFranky(isCoulomb)
nDOF = 6;

% orientation
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
XYZ = eye(3);
mXZY = [-ex ez ey];
XmZY = [ex -ez ey];
mXmYZ = [-ex -ey ez];
ZYmX = [ez ey -ex];
XZmY = [ex ez -ey];
mYmXmZ = [-ey -ex -ez];

% LINK FRAMES
baseFrame = eye(4);
baseFrame(1:3,1:3) = XYZ; baseFrame(1:3,4) = 0.001 * [0;-8;0];
linkFrame = eye(4);
initialLinkFrames = cell(nDOF,1);
linkFrame(1:3,1:3) = XYZ; linkFrame(1:3,4) = 0.001 * [0;87;0];
initialLinkFrames{1} = linkFrame;
linkFrame(1:3,1:3) = mXZY; linkFrame(1:3,4) = 0.001 * [0;152;63];
initialLinkFrames{2} = linkFrame;
linkFrame(1:3,1:3) = XmZY; linkFrame(1:3,4) = 0.001 * [0;396;29];
initialLinkFrames{3} = linkFrame;
linkFrame(1:3,1:3) = XYZ; linkFrame(1:3,4) = 0.001 * [0;641;-29];
initialLinkFrames{4} = linkFrame;
linkFrame(1:3,1:3) = ZYmX; linkFrame(1:3,4) = 0.001 * [0;627.25;131];
initialLinkFrames{5} = linkFrame;
linkFrame(1:3,1:3) = XZmY; linkFrame(1:3,4) = 0.001 * [0;694;155.50];
initialLinkFrames{6} = linkFrame;

% MOTOR FRAMES
motorFrame = eye(4);
initialMotorFrames = cell(nDOF,1);
motorFrame(1:3,1:3) = XmZY; motorFrame(1:3,4) = 0.001 * [0; 81; 0];
initialMotorFrames{1} = motorFrame;
motorFrame(1:3,1:3) = XYZ; motorFrame(1:3,4) = 0.001 * [0; 152; 57];
initialMotorFrames{2} = motorFrame;
motorFrame(1:3,1:3) = mYmXmZ; motorFrame(1:3,4) = 0.001 * [0; 396; 35];
initialMotorFrames{3} = motorFrame;
motorFrame(1:3,1:3) = XYZ; motorFrame(1:3,4) = 0.001 * [0; 609; -29];
initialMotorFrames{4} = motorFrame;
motorFrame(1:3,1:3) = mXZY; motorFrame(1:3,4) = 0.001 * [0; 602.75; 131];
initialMotorFrames{5} = motorFrame;
motorFrame(1:3,1:3) = mXmYZ; motorFrame(1:3,4) = 0.001 * [0; 694; 131];
initialMotorFrames{6} = motorFrame;
nMotor = length(initialMotorFrames);

% JOINT SCREW w.r.t. motor frame
temp = cell2mat(initialMotorFrames);
qtemp = reshape(temp(:,4),4,[]);
wtemp = reshape(temp(:,3),4,[]);
q_screw = qtemp(1:3,:);
w_screw = wtemp(1:3,:);
initialScrew = zeros(6,nDOF);
for i=1:nDOF
    % screw
    initialScrew(1:3,i) = w_screw(:,i);
    initialScrew(4:6,i) = skew(q_screw(:,i)) * w_screw(:,i);
end

% MOTOR INERTIA
phi_PH42 = getInertiaMotorPH42();
phi_PH54 = getInertiaMotorPH54();
phi_XH540 = getInertiaMotorXH540();
% All motor inertia w.r.t. reference frame
Phi_motors = [phi_PH42 phi_PH54 phi_PH42 phi_XH540 phi_XH540 phi_XH540];
% coordiante transform: motor coordinates to reference frame
motorPhi = cell(nMotor,1);
for i=1:nMotor
    motorPhi{i} = getInertiaTransformRegressor(initialMotorFrames{i}) * Phi_motors(:,i);
end

% MOTOR LIMITS
% velocity (rad/sec)
rpm = 2*pi/60;
maxVel_PH42 = rpm * 2920 * 0.01;
maxVel_PH54 = rpm * 2900 * 0.01; % range: 0 ~ 1,023
maxVel_XH54 = rpm * 128 * 0.229; % range: 0 ~ 1,023
motorVelLimit = [maxVel_PH42 maxVel_PH54 maxVel_PH42 maxVel_XH54 maxVel_XH54 maxVel_XH54];
% current (Ampere) 
millierAmpere = 0.001;
maxCurr_PH42 = 4500 * millierAmpere;
maxCurr_PH54 = 22740 * millierAmpere;
maxCurr_XH54 = 2047 * 2.69*millierAmpere;
motorCurrentLimit = [maxCurr_PH42 maxCurr_PH54 maxCurr_PH42 maxCurr_XH54 maxCurr_XH54 maxCurr_XH54];

% Motor current gain ( joint torque = gain * current)
gain_PH42 = 4.69; % Nm / Ampere
gain_PH54 = 5.28; % Nm / Ampere
gain_XH54 = 4.17; % Nm / Ampere
currentGain = [gain_PH42 gain_PH54 gain_PH42 gain_XH54 gain_XH54 gain_XH54];

% End-effector frame
M_initial = initialLinkFrames{end};
M_initial(1:3,1:3) = eye(3);
M_initial(3,4) = M_initial(3,4) + 0.07;

% struct
initialRobot = struct;
initialRobot.baseFrame = baseFrame;
initialRobot.linkFrames= initialLinkFrames;
initialRobot.motorFrames = initialMotorFrames;
initialRobot.screw = initialScrew;
initialRobot.motorPhi = motorPhi;
initialRobot.coulomb = diag(zeros(6,1));

if isCoulomb
    initialRobot.coulomb = diag([0.505522078192190   1.654455354978156   0.654146578009750   0.263549934816684   0.026795649657167   0.025090728181183]);
    initialRobot.viscous = diag([-0.073135860905893   2.619584501483819   0.088178887265450   0.192854684880986   0.128443069189294   0.117764920830548]);
    initialRobot.rotorinertia = diag([0.092225988589875   0.022611181425337  -0.033497602750333  -0.010869828291476   0.022104793279946   0.026610457227956]);
else
    initialRobot.coulomb = diag(zeros(1,nDOF));
    initialRobot.viscous = diag([0.289091153640209   5.420179095964599   0.953562870084288   0.489511755071330   0.212032026197624   0.130307445371531]);
    initialRobot.rotorinertia = diag([0.136913464387272   0.170966917431287   0.130824945335699   0.028212197285790   0.029043364168807   0.039129572056044]);
end

initialRobot.k = ones(1,nDOF);
initialRobot.motorPosLimit = [0.8*pi, 0.5*pi, 0.8*pi, 0.6*pi, 2*pi, 2*pi];
initialRobot.motorVelLimit = motorVelLimit;
initialRobot.motorCurrentLimit = motorCurrentLimit;
initialRobot.currentGain = currentGain; % joint torque 
initialRobot.payload = 0;
initialRobot.M = M_initial;

end
