function meshGroup = getMeshGroupFranky(linkMeshNames, motorMeshNames, initialRobot)
nLink = length(linkMeshNames);
density_ABS = 1.05 * 1e-3 / 1e-6 * 0.7971; % measured

if length(linkMeshNames) ~= length(motorMeshNames) || length(linkMeshNames) ~= length(initialRobot.linkFrames)
    linkMeshNames
    motorMeshNames
    initialRobot.linkFrames
    error('Number of the links must be the same as the motors.')
end

% VERTEX SEGMENT
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
expandingDirection = zeros(3,nLink); % 1=x, 2=y, 3=z
c_begin = zeros(1,nLink);
c_end= zeros(1,nLink);
% link 1
expandingDirection(:,1) = ey;
c_begin(1) = 1e-3 * 95;
c_end(1) = 1e-3 * 117;
% link 2
expandingDirection(:,2) = ey;
c_begin(2) = 1e-3 * 180;
c_end(2) = 1e-3 * 360;
% link 3
expandingDirection(:,3) = ey;
c_begin(3) = 1e-3 * 418;
c_end(3) = 1e-3 * 544;
% link 4
expandingDirection(:,4) = ez;
c_begin(4) = 1e-3 * (-2);
c_end(4) = 1e-3 * 69;
% link 5
expandingDirection(:,5) = ey;
c_begin(5) = 1e-3 * 634;
c_end(5) = 1e-3 * 674;
% link 6
expandingDirection(:,6) = ez;
c_begin(6) = 1e-3 * 163;
c_end(6) = 1e-3 * 200;
% adjustable link length
initialLength = c_end - c_begin;
if any(initialLength<0)
    error(['c_begin must be less than c_end! c_end<c_begin=' num2str(c_end<c_begin)])
end

% VERTEX GROUP
meshGroup = {};
RotX = expm(skew([1;0;0])*pi/2);
for i=1:nLink
    % MOTOR
    [face, vert] = plyread(motorMeshNames{i},'tri');
    vert = vert * 1e-3; % millimeter to meter
    if i <= 3
        vert = vert * RotX'; % Pro-series motor axis = z-axis. (but y-axis in CAD)
    end
    transformed4vec = initialRobot.motorFrames{i} * [vert';ones(1,size(vert,1))];
    vert = transformed4vec(1:3,:)';
    mesh = struct;
    mesh.vertices = vert;
    mesh.faces = face;
    mesh.isLink = false;
    mesh.frame = initialRobot.motorFrames{i};
    mesh.parent = length(meshGroup); % parent mesh
    mesh.motorNo = i;
    meshGroup{length(meshGroup)+1} = mesh;
    
    % LINK
    [face, vert] = plyread(linkMeshNames{i},'tri');
    vert = vert * 1e-3; % millimeter to meter
    transformed4vec = initialRobot.linkFrames{i} * [vert';ones(1,size(vert,1))];
    vert = transformed4vec(1:3,:)';
    % segment
    c = vert * expandingDirection(:,i);
    id_before = c <= c_begin(i);
    id_after = c_end(i) <= c;
    id_adjust = ~id_before & ~id_after;    
    mesh = struct;
    mesh.vertices = vert;
    mesh.faces = face;
    mesh.isLink = true;
    mesh.parent = length(meshGroup); % parent mesh
    mesh.frame = initialRobot.linkFrames{i};
    mesh.linkNo = i;
    mesh.id_before = id_before;
    mesh.id_adjust = id_adjust;
    mesh.id_after = id_after;
    mesh.expandDirection = expandingDirection(:,i);
    mesh.initialLength = initialLength(i);
    mesh.A_vertex = expandingDirection(:,i) * expandingDirection(:,i)';
    mesh.b_vertex = c_begin(i) * expandingDirection(:,i);
    mesh.density = density_ABS;
    meshGroup{length(meshGroup)+1} = mesh;
end

end