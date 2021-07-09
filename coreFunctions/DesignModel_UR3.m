function [A_screw, M_screw, Phi, dEtaA_dRho, dEtaM_dRho, dPhi_dRho, meshGroup, TransformSE3] = DesignModel_UR3(meshGroup_initial, designParameters, initialRobot, getGradient)

% % INPUTS
% meshGroup_initial: cell of the initial trimeshes
% designParameters: the design parameters (rho)
% initialRobot: robot initial configurations
% getGradient: boolean to request the derivatives

% % OUTPUTS
% A_screw: joint screws (6 x n) 
% M_screw: log of the end-effector pose (6 x 1)
% Phi: Concatenated vector of inertia (10n x 1)
% dEtaA_dRho: screw derivative w.r.t design parameter
% dEtaM_dRho: end-effector pose derivative w.r.t design parameter
% dPhi_dRho: inertia derivative w.r.t design parameter
% meshGroup: cell of the deformed trimeshes
% TransformSE3: Rigid body motions induced by each link deformation

nDesignParam = length(designParameters);
nLink = size(initialRobot.screw,2);
if any(size(designParameters)~=[nLink,1]) || length(meshGroup_initial) ~= 2*nLink
    error(['DesignModel_UR3(): size(designParameters)=' num2str(size(designParameters))])
end
if nargin < 4
    getGradient = true;
end
    
% Kinematic parameters
meshGroup = meshGroup_initial;
[A_screw, M_screw, dEtaA_dRho, dEtaM_dRho, TransformSE3, spaceJacobian_rho] = DesignModel_UR3Kinematics(meshGroup_initial, designParameters, initialRobot, getGradient);

% Vertex transformation
exp_rho = exp(designParameters);
count_transform = 0;
LinkNoToCount_transform = zeros(nLink,1);
LinkNoToMeshNo = zeros(nLink,1);
diffVertexByDesignParam = cell(nLink,1); % each cell = (3*nVertex) x nDesignParam
for i=1:length(meshGroup_initial)
    mesh = meshGroup_initial{i};
    if mesh.isLink
        count_transform = count_transform + 1;
        LinkNoToMeshNo(mesh.linkNo) = i;
        LinkNoToCount_transform(mesh.linkNo) = count_transform;
        nVert = size(mesh.vertices,1);
        if count_transform == 1
            SE3 = eye(4);
        else
            SE3 = TransformSE3{count_transform-1};
        end
        R = SE3(1:3,1:3);
        p = SE3(1:3,4);
        % before
        mesh.vertices(mesh.id_before,:) = (R * mesh.vertices(mesh.id_before,:)' + p)';
        mesh.frame = SE3 * mesh.frame;
        % adjust
        R_Av_b = R * (mesh.A_vertex * mesh.vertices(mesh.id_adjust,:)' - mesh.b_vertex);
        mesh.vertices(mesh.id_adjust,:) = (R * mesh.vertices(mesh.id_adjust,:)' + p + (exp_rho(count_transform)-1)*R_Av_b)';
        % after
        SE3 = TransformSE3{count_transform};
        R = SE3(1:3,1:3);
        p = SE3(1:3,4);
        mesh.vertices(mesh.id_after,:) = (R * mesh.vertices(mesh.id_after,:)' + p)';
        % differential TransformSE3 = Space jacobian w.r.t. the design parameter
        % differential vertex
        if getGradient
            diffVertex = zeros(3*nVert, nDesignParam);
            total_vertex = reshape(1:(3*nVert),3,[])';
            total_adjust = reshape(total_vertex(mesh.id_adjust,:)',[],1);
            total_after = reshape(total_vertex(mesh.id_after,:)',[],1);
            % before + adjust + after
            if count_transform > 1
                for j=1:count_transform-1
                    w = spaceJacobian_rho(1:3,j);
                    v = spaceJacobian_rho(4:6,j);
                    diffVertex(:,j) = reshape(v + skew(w) * mesh.vertices',[],1);
                end
            end
            % adjust
            diffVertex(total_adjust,count_transform)= diffVertex(total_adjust,count_transform) + reshape(exp_rho(count_transform) * R_Av_b,[],1);
            w = spaceJacobian_rho(1:3,count_transform);
            v = spaceJacobian_rho(4:6,count_transform);
            % after
            diffVertex(total_after,count_transform)= diffVertex(total_after,count_transform) + reshape(v + skew(w) * mesh.vertices(mesh.id_after,:)',[],1);
            diffVertexByDesignParam{mesh.linkNo} = diffVertex;
        end
    else % motor
        if count_transform > 0
            SE3 = TransformSE3{count_transform};
            R = SE3(1:3,1:3);
            p = SE3(1:3,4);
            mesh.vertices = (R * mesh.vertices' + p)';
            mesh.frame = SE3 * mesh.frame;
        end
    end
    meshGroup{i} = mesh;
end
if count_transform ~= nDesignParam
    error(['count_transform must be length(designParameters) at the end: count_transform=' num2str(count_transform) ', length(designParameters)=' num2str(nDesignParam)])
end

% INERTIA (links by mesh, motors by manufacturer-provided values)
Phi = zeros(10*nLink,1);
dPhi_dRho = zeros(10*nLink, nDesignParam);
% inertia
for i=1:nLink
    id_phi = 10*(i-1)+1:10*i;
    mesh = meshGroup{LinkNoToMeshNo(i)};
    if i < nLink
        SE3 = TransformSE3{LinkNoToCount_transform(i)}; 
        InvAdj = largeAdjoint(inverseSE3(SE3));
        phi_motor = GToPhi_latest(InvAdj' * PhiToG_latest(initialRobot.motorPhi{i+1}) * InvAdj);
    else % i == nLink
        phi_motor = zeros(10,1); % payload
    end
    if getGradient
        [phi_link, dPhi_link] = MeshToInertia(mesh.vertices, mesh.faces, mesh.density);
    else
        phi_link = MeshToInertia(mesh.vertices, mesh.faces, mesh.density);
    end
    
    Phi(id_phi) = phi_link + phi_motor;
    
    % differential inertia
    if getGradient
        for j = 1:LinkNoToCount_transform(i)
            dPhi_dRho(id_phi, j) = translateInertia(phi_motor, spaceJacobian_rho(:,j)); % differential motor inertia by translation
        end
        dPhi_dRho(id_phi,:) = dPhi_dRho(id_phi,:) + dPhi_link * diffVertexByDesignParam{i}; % differential link inertia by differential vertex
    end
end

end
