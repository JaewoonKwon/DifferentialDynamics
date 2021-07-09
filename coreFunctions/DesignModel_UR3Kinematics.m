function [A_screw, M_screw, dEtaA_dRho, dEtaM_dRho, TransformSE3, spaceJacobian_rho] = DesignModel_UR3Kinematics(meshGroup_initial, designParameters, initialRobot, getGradient)

% % INPUTS
% meshGroup_initial: cell of the initial trimeshes
% designParameters: the design parameters (rho)
% initialRobot: robot initial configurations
% getGradient: boolean to request the derivatives

% % OUTPUTS
% A_screw: joint screws (6 x n) 
% M_screw: log of the end-effector pose (6 x 1)
% dEtaA_dRho: screw derivative w.r.t design parameter
% dEtaM_dRho: end-effector pose derivative w.r.t design parameter
% TransformSE3: Rigid body motions induced by each link deformation
% spaceJacobian_rho: body Jacobian the link deformations

nDesignParam = length(designParameters);
nLink = size(initialRobot.screw,2);
if any(size(designParameters)~=[nLink,1]) || length(meshGroup_initial) ~= 2*nLink
    error(['DesignModel_UR3Kinematics(): size(designParameters)=' num2str(size(designParameters))])
end
if nargin < 4
    getGradient = true;
end

% Initial end-effector pose
M_initial = initialRobot.M;

% Vertex translation by design parameter
TransformSE3 = cell(nDesignParam,1);
spaceJacobian_rho = zeros(6,nDesignParam);
exp_rho = exp(designParameters);
rhoSE3 = eye(4);
LinkNoTocount_transform = zeros(nLink,1);
count_transform = 0;
for i=1:length(meshGroup_initial)
    mesh = meshGroup_initial{i};
    if mesh.isLink
        count_transform = count_transform + 1;
        LinkNoTocount_transform(mesh.linkNo) = count_transform;
        if getGradient
            v_rho = mesh.expandDirection * mesh.initialLength * exp_rho(count_transform);
        end
        if count_transform == 1
            rhoSE3(1:3,4) = mesh.expandDirection * mesh.initialLength * (exp_rho(count_transform)-1);
            TransformSE3{count_transform} = rhoSE3;
            if getGradient
                spaceJacobian_rho(:,1) = [zeros(3,1);v_rho];
            end
        else
            rhoSE3(1:3,4) = mesh.expandDirection * mesh.initialLength * (exp_rho(count_transform)-1);
            TransformSE3{count_transform} = TransformSE3{count_transform-1} * rhoSE3;
            if getGradient
                spaceJacobian_rho(:,count_transform) = largeAdjoint(TransformSE3{count_transform-1}) * [zeros(3,1);v_rho];
            end
        end
    end
end

% SCREW
A_screw = zeros(6,nLink);
dEtaA_dRho = zeros(6*nLink, nDesignParam);
for i=1:nLink
    id_screw = 6*(i-1)+1:6*i;
    if i == 1
        Adj = eye(6);
    else
        SE3 = TransformSE3{LinkNoTocount_transform(i-1)};
        Adj = largeAdjoint(SE3);
    end
    % screw
    screw = initialRobot.screw(:,i);
    A_screw(:,i) = Adj * screw;    
    % differential screw
    if i>1 && getGradient
        dEtaA_dRho(id_screw, 1:LinkNoTocount_transform(i-1)) = spaceJacobian_rho(:,1:LinkNoTocount_transform(i-1));
    end
end
% End-effector pose of the new design
M_screw = ToVector(logm(TransformSE3{end} * M_initial));
% if getGradient
    dEtaM_dRho = spaceJacobian_rho;
% end
end
