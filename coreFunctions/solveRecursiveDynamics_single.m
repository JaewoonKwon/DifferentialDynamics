function state = solveRecursiveDynamics_single(pos, vel, acc, Phi, A_screw, initialLinkFrames, Vdot0, F_ext)

%  (Abbreviation: "ExpBodyFr" = each expressed in the respective link body frame.)

% % INPUTS
% pos: joint angle (1 x n)
% vel: joint velocity (1 x n)
% acc: joint acceleration (1 x n)
% Phi: inertia (ExpBodyFr, 10n x 1)
% A_screw: joint twists (ExpBodyFr, 6 x n)
% initialLinkFrames: initial link frames (at the zero configuration)
% Vdot0: acceleration of the base (ExpBodyFr, only gravity if fixed-base)
% F_ext: contact force wrench exerted on the end-effector (expr. in the last link frame)

% % OUTPUTS
% "state" containing the followings:
% jointPos: joint angle (same as input 'pos')
% jointVel: joint velocity (same as input 'vel')
% jointAcc: joint acceleration (same as input 'acc')
% linkFrames: initial link frames (same as input 'initialLinkFrames')
% V: link body velocities (ExpBodyFr, 6 x n)
% Vdot: link body accelerations (ExpBodyFr, 6 x n)
% F: link body force wrenches (ExpBodyFr, 6 x n)
%      ((i-th col.) F_i denotes the force exerted on link i by link i-1) 
% W: regressor of the link forces (6n x 10n)
% Y: regressor of the joint torques (6 x 10n)
% jointTorque: joint torque (1 x n)

% exception
nJoint = size(A_screw,2);
if any(size(pos) ~= [1, nJoint]) || any(size(vel) ~= [1, nJoint]) || any(size(acc) ~= [1, nJoint])
    size(pos)
    size(vel)
    size(acc)
    error('Error 1')
end
if any(size(A_screw) ~= [6, nJoint]) || length(initialLinkFrames) ~= nJoint
    size(A_screw)
    size(initialLinkFrames)
    error('Error 2')
end
if any(size(Vdot0) ~= [6, 1]) ||any(size(F_ext) ~= [6, 1])
    size(Vdot0)
    size(F_ext)
    error('Error 3')
end
% variables
linkFrames = cell(size(initialLinkFrames));
V = zeros(6,nJoint);
Vdot = zeros(6,nJoint);
F = zeros(6,nJoint);
jointTorque = zeros(1, nJoint);

% RECURSIVE DYNAMICS
inertia = reshape(Phi,10,[]);
G = zeros(6,6,nJoint);
for i=1:nJoint
    G(:,:,i) = PhiToG_latest(inertia(:,i));
end
% Forward: T, V, Vdot, P_cell, Q_cell
V_parent = zeros(6,1);
Vdot_parent = Vdot0;
for i=1:nJoint
    A_i = A_screw(:,i);
    % V, Vdot, T
    T_i = initialLinkFrames{i} * expm(ToMatrix(A_i)*pos(i));
    AdjInvT = largeAdjoint(inverseSE3(T_i));
    V_i = A_i * vel(i) + AdjInvT * V_parent;
    Vdot_i = A_i * acc(i) + AdjInvT * Vdot_parent + smallAdjoint(V_i) * A_i * vel(i);
    % log
    V_parent = V_i;
    Vdot_parent = Vdot_i;
    linkFrames{i} = T_i;
    V(:,i) = V_i;
    Vdot(:,i) = Vdot_i;
end
% Backward: F, torque, R_cell
F_child = zeros(6,1);
W = zeros(6 * nJoint, 10 * nJoint);
Y = zeros(nJoint, 10 * nJoint);
for i=nJoint:-1:1
    A_i = A_screw(:,i);
    G_i = G(:,:,i);
    V_i = V(:,i);
    Vdot_i = Vdot(:,i);
    if i<nJoint
        T_ip1 = linkFrames{i+1};
    else
        T_ip1 = eye(4);
    end
    AdjInvT_ip1 = largeAdjoint(inverseSE3(T_ip1));
    adV_i = smallAdjoint(V_i);
    % Y, W (F = W * Phi)
    if i==nJoint
        W_diagonal = V_regressor_latest(Vdot_i) - adV_i' * V_regressor_latest(V_i);
        W(6*(i-1)+1:6*i,10*(i-1)+1:end) = W_diagonal;
    else
        W_diagonal = V_regressor_latest(Vdot_i) - adV_i' * V_regressor_latest(V_i);
        W_offdiag  = transpose(AdjInvT_ip1) * W(6*i+1:6*(i+1),10*i+1:end);
        W(6*(i-1)+1:6*i,10*(i-1)+1:end) = [W_diagonal W_offdiag];
    end
    Y(i,:) = A_i' * W(6*(i-1)+1:6*i,:);
    % F, torque
    if i==nJoint
        F_i = G_i * Vdot_i - adV_i' * G_i * V_i - F_ext;
    else
        F_i = G_i * Vdot_i - adV_i' * G_i * V_i + AdjInvT_ip1' * F_child;
    end
    tau_i = A_i' * F_i;
    % log
    F_child = F_i;
    F(:,i) = F_i;
    jointTorque(i) = tau_i;
end

% output
state = struct;
state.jointPos = pos;
state.jointVel = vel;
state.jointAcc = acc;
state.linkFrames = linkFrames;
state.V = V;
state.Vdot = Vdot;
state.F = F;
state.W = W;
state.Y = Y;
state.jointTorque = jointTorque;
end