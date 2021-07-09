% computation time
clear all
close all
clc

addpath('mesh');
addpath('basicFunctions');
addpath('coreFunctions');

nCompute = 100;
nDisp = round(nCompute/10);
disp('--------- Function computation time ---------')
useSimplifiedMesh = true;

%% common variables
initialRobot = getRobot(true);
nJoint = size(initialRobot.screw,2);
Vdot0 = [0 0 0 0 9.81 0]';
initialLinkFrames = mat2cell(repmat(eye(4),nJoint,1),ones(nJoint,1)*4, 4);
F_ext = zeros(1,6);
if useSimplifiedMesh
    linkMeshNames = {'link1.ply', 'link2.ply', 'link3.ply', 'link4.ply', 'link5.ply', 'link6.ply'};
else
    linkMeshNames = {'link1_sw.ply', 'link2_sw.ply', 'link3_sw.ply', 'link4_sw.ply', 'link5_sw.ply', 'link6_sw.ply'};
end

motorMeshNames = {'H42P-020-S300-R.ply', 'H54P-200_M54P-060.ply', 'H42P-020-S300-R.ply', 'XH-540_idle.ply', 'XH-540_idle.ply', 'XH-540_idle.ply'};
meshGroup_initial = getMeshGroupFranky(linkMeshNames, motorMeshNames, initialRobot);
% initial robot design
designParameters = rand(nJoint,1);
[A_screw, M_screw, Phi, ~, ~, ~, ~] = DesignModel_UR3(meshGroup_initial, designParameters, initialRobot);

thetaPos = rand(nCompute,nJoint);
thetaVel = rand(nCompute,nJoint);
thetaAcc = rand(nCompute,nJoint);

nVertices = zeros(1, nJoint);
for i=1:length(meshGroup_initial)
    mesh = meshGroup_initial{i};
    if mesh.isLink
        nVert = size(mesh.vertices,1);
        nFace = size(mesh.faces,1);
        nVertices(mesh.linkNo) = nVert;
        disp([num2str(mesh.linkNo) 'th link: # of vertices = ' num2str(nVert) ',  # of faces = ' num2str(nFace)])
    end
end

%% Recursive Dynamics
timeRecurDyn = zeros(nCompute,1);
for i=1:nCompute
    tic
    state = solveRecursiveDynamics_single(thetaPos(i,:), thetaVel(i,:), thetaAcc(i,:), Phi, A_screw, initialLinkFrames, Vdot0, F_ext);
    timeRecurDyn(i) = toc;
    if mod(i,nDisp) == 0
        disp(['[Recursive Dynamics] count/max = ' num2str(i) '/' num2str(nCompute)])
    end
end

% Differnetial Dynamics
timeDiffDyn = zeros(nCompute,1);
for i=1:nCompute
    tic
    state = solveDifferentialDynamics_single(thetaPos(i,:), thetaVel(i,:), thetaAcc(i,:), Phi, A_screw, initialLinkFrames, Vdot0, F_ext);
    timeDiffDyn(i) = toc;
    if mod(i,nDisp) == 0
        disp(['[Differential Dynamics] count/max = ' num2str(i) '/' num2str(nCompute)])
    end
end

%% Design model
getGradient = false;
timeDesignModel = zeros(nCompute,1);
for i=1:nCompute
    tic
    [A_screw, M_screw, Phi, dEtaA_dRho, dEtaM_dRho, dPhi_dRho, meshGroup, TransformSE3] = DesignModel_UR3(meshGroup_initial, designParameters, initialRobot, getGradient);
    timeDesignModel(i) = toc;
    if mod(i,nDisp) == 0
        disp(['[Design Model] count/max = ' num2str(i) '/' num2str(nCompute)])
    end
end

% Differential Design Model
getGradient = true;
timeDiffDesignModel = zeros(nCompute,1);
for i=1:nCompute
    tic
    [A_screw, M_screw, Phi, dEtaA_dRho, dEtaM_dRho, dPhi_dRho, meshGroup, TransformSE3] = DesignModel_UR3(meshGroup_initial, designParameters, initialRobot, getGradient);
    timeDiffDesignModel(i) = toc;
    if mod(i,nDisp) == 0
        disp(['[Differential Design Model] count/max = ' num2str(i) '/' num2str(nCompute)])
    end
end

%% Result

disp('---------------- SUMMARY ----------------')
disp(['Recursive Dynamics: mean = ' num2str(mean(timeRecurDyn)) ' (' num2str(1/mean(timeRecurDyn)) ' Hz), std = ' num2str(std(timeRecurDyn))])
disp(['Differential Dynamics: mean = ' num2str(mean(timeDiffDyn)) ' (' num2str(1/mean(timeDiffDyn)) ' Hz), std = ' num2str(std(timeDiffDyn))])
disp(['Design Model: mean = ' num2str(mean(timeDesignModel)) ' (' num2str(1/mean(timeDesignModel)) ' Hz), std = ' num2str(std(timeDesignModel))])
disp(['Differential Design Model: mean = ' num2str(mean(timeDiffDesignModel)) ' (' num2str(1/mean(timeDiffDesignModel)) ' Hz), std = ' num2str(std(timeDiffDesignModel))])


