function state = solveDifferentialDynamics_single(pos, vel, acc, Phi, A_screw, initialLinkFrames, Vdot0, F_ext)

% % INPUTS
% pos: joint angle (1 x n)
% vel: joint velocity (1 x n)
% acc: joint acceleration (1 x n)
% Phi: inertia (10n x 1)
% A_screw: joint twists (6 x n)
% initialLinkFrames: initial link frames (at the zero configuration)
% Vdot0: acceleration of the base (only gravity if fixed-base)
% F_ext: contact force wrench exerted on the end-effector

% % OUTPUTS
% "state" containing the followings:
% jointPos: joint angle (same as input 'pos')
% jointVel: joint velocity (same as input 'vel')
% jointAcc: joint acceleration (same as input 'acc')
% linkFrames: initial link frames (same as input 'initialLinkFrames')
% V: link body velocities (6 x n)
% Vdot: link body accelerations (6 x n)
% F: link body force wrenches (6 x n)
% W: regressor of the link forces (6n x 10n)
% Y: regressor of the joint torques (6 x 10n)
% jointTorque: joint torque (1 x n)
% P_cell: the derivatives of V 
% Q_cell: the derivatives of Vdot
% R_cell: the derivatives of F
% S: the derivatives of the joint torque

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
if any(size(Vdot0) ~= [6, 1]) ||any(size(F_ext) ~= [1, 6])
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
P_cell = cell(nJoint,1);
Q_cell = cell(nJoint,1);
for i=1:nJoint
    A_i = A_screw(:,i);
    % V, Vdot, T
    T_i = initialLinkFrames{i} * expm(ToMatrix(A_i)*pos(i));
    AdjInvT = largeAdjoint(inverseSE3(T_i));
    V_i = A_i * vel(i) + AdjInvT * V_parent;
    Vdot_i = A_i * acc(i) + AdjInvT * Vdot_parent + smallAdjoint(V_i) * A_i * vel(i);
    % P, Q
    Pr = getPr(i, nJoint, A_i, V_parent, T_i, vel);
    Qr = getQr(i, nJoint, A_i, Vdot_parent, V_i, T_i, vel, acc);
    if i==1
        P_cell{i} = Pr;
        Q_cell{i} = Qr - smallAdjoint(A_i*vel(i))*P_cell{i};
    else
        P_cell{i} = AdjInvT * P_cell{i-1} + Pr;
        Q_cell{i} = AdjInvT * Q_cell{i-1} + Qr - smallAdjoint(A_i*vel(i))*P_cell{i};
    end
    % log
    V_parent = V_i;
    Vdot_parent = Vdot_i;
    linkFrames{i} = T_i;
    V(:,i) = V_i;
    Vdot(:,i) = Vdot_i;
end
% Backward: F, torque, R_cell
F_child = F_ext';
R_cell = cell(nJoint,1);
S = zeros(nJoint, 9*nJoint);
W = zeros(6 * nJoint, 10 * nJoint);
Y = zeros(nJoint, 10 * nJoint);
for i=nJoint:-1:1
    A_i = A_screw(:,i);
    if i<nJoint
        A_ip1 = A_screw(:,i+1);
    else
        A_ip1 = zeros(6,1);
    end
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
    F_i = AdjInvT_ip1' * F_child + G_i * Vdot_i - adV_i' * G_i * V_i;
    tau_i = A_i' * F_i;
    % R
    Rr = getRr(i, nJoint, T_ip1, F_child, A_ip1);
    if i==nJoint
        R_cell{i} = Rr + (-adV_i'*G_i - smallAd_star(G_i*V_i))*P_cell{i} + G_i*Q_cell{i};
    else
        R_cell{i} = AdjInvT_ip1' * R_cell{i+1} + Rr + (-adV_i'*G_i - smallAd_star(G_i*V_i))*P_cell{i} + G_i*Q_cell{i};
    end
    % S
    S(i,:) = A_i'*R_cell{i};
    S(i,6*(i-1)+1:6*i) = S(i,6*(i-1)+1:6*i) - F_i' * smallAdjoint(A_i);
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
state.P_cell = P_cell;
state.Q_cell = Q_cell;
state.R_cell = R_cell;
state.S = S;
end

function Pr = getPr(i, nJoint, A_i, V_parent, T_i, vel)
Pr = zeros(6,9*nJoint);
AdjT = largeAdjoint(T_i);
AdjInvT = largeAdjoint(inverseSE3(T_i));
adV_parent = smallAdjoint(V_parent);
Pr(:,1+6*(i-1):6*i) = AdjInvT * adV_parent * (eye(6)-AdjT) - smallAdjoint(A_i*vel(i));
Pr(:,i+6*nJoint) = AdjInvT * adV_parent * A_i;
Pr(:,i+7*nJoint) = A_i;
end

function Qr = getQr(i, nJoint, A_i, Vdot_parent, V_i, T_i, vel, acc)
Qr = zeros(6,9*nJoint);
AdjT = largeAdjoint(T_i);
AdjInvT = largeAdjoint(inverseSE3(T_i));
adVdot_parent = smallAdjoint(Vdot_parent);
adA_i = smallAdjoint(A_i);
adV_i = smallAdjoint(V_i);
Qr(:,1+6*(i-1):6*i) = AdjInvT * adVdot_parent * (eye(6)-AdjT) - (adV_i*vel(i) + acc(i)*eye(6)) * adA_i;
Qr(:,i+6*nJoint) = AdjInvT * adVdot_parent * A_i;
Qr(:,i+7*nJoint)  = adV_i * A_i;
Qr(:,i+8*nJoint)  = A_i;
end

function Rr = getRr(i, nJoint, T_ip1, F_child, A_ip1)
Rr = zeros(6,9*nJoint);
AdjInvT_ip1 = largeAdjoint(inverseSE3(T_ip1));
TransAdjInvT_ip1 = AdjInvT_ip1';
adstarF_child = smallAd_star(F_child);
if i<nJoint
    Rr(:,1+6*i:6*(i+1)) = - TransAdjInvT_ip1 * adstarF_child * (AdjInvT_ip1-eye(6));
    Rr(:,i+1+6*nJoint) = - TransAdjInvT_ip1 * adstarF_child * A_ip1;
end
end

