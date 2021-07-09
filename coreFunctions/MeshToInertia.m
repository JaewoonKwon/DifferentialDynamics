function [phi, dPhidVertex] = MeshToInertia(vertices, face, density)

% % INPUTS
% vertices: the mesh vertices (nVert x 3)
% face: the mesh faces (nFace x 3)
% density: density of the link

% % OUTPUTS
% phi = mass, CoM, inertia (in the following order: xx, yy, zz, xy, xz, yz) (10 x 1)
% dPhidVertex: phi derivative w.r.t vertices (10 x 3nVert)

nVert = size(vertices,1);
nFace = size(face,1);
if size(vertices,2) ~= 3 || size(face,2) ~= 3
    error(['MeshToInertia(): Wrong input. size(xyz)=' num2str(size(vertices)) ', size(normal)=' num2str(size(face))])
end

m_multiplied18 = 0; % mass
h_multiplied24 = zeros(3,1); % mass x CoM
rotInert_diag = zeros(3,1); % diagonal of rotational inertia
rotInert_off = zeros(3,1); % off-diagonal of rotational inertia
if nargout >= 2
    dm_multiplied18 = zeros(1,3*nVert);
    dh_multiplied24 = zeros(3,3*nVert);
    drotInert_diag = zeros(3,3*nVert);
    drotInert_off = zeros(3,3*nVert);
end
for i=1:nFace
    % vertices
    id_vert = face(i,:);
    vert = vertices(id_vert,:);
    vi = vert(1,:);
    vj = vert(2,:);
    vk = vert(3,:);
    e1 = vj - vi;
    e2 = vk - vi;
    % Area (normal) vector
    normal = e2 * (-skew(e1)); % transpose(skew(e1)*e2')
    v_ijk = vi + vj + vk;
    if nargout == 1
        g = g_vertex(vi, vj, vk);
    elseif nargout >= 2
        [g, dg] = g_vertex(vi, vj, vk);
    end
    % mass
    m_multiplied18 = m_multiplied18 + normal * v_ijk';
    % center of mass (multiplied by mass)
    h_multiplied24 = h_multiplied24 + (g .* normal');
    % rotational inertia
    det_J = det(vert);
    rot_g_diag = [g(2) + g(3); g(1) + g(3); g(1) + g(2)];
    % diagonal
    rotInert_diag = rotInert_diag + det_J * rot_g_diag;
    % off-diagonal
    r_ijk = r(vi)*vi' + r(vi)*vj' + r(vi)*vk' + r(vj)*vj' + r(vj)*vk' + r(vk)*vk';
    rotInert_off = rotInert_off - det_J * r_ijk;
    if nargout >= 2
        % vertice index
        id_i = 3*(id_vert(1)-1)+1:3*id_vert(1);
        id_j = 3*(id_vert(2)-1)+1:3*id_vert(2);
        id_k = 3*(id_vert(3)-1)+1:3*id_vert(3);
        % mass
        dnormal_i = skew(vk-vj);
        dnormal_j = skew(vi-vk);
        dnormal_k = skew(vj-vi);
        dm_multiplied18(:,id_i) = dm_multiplied18(:,id_i) + v_ijk * dnormal_i + normal;
        dm_multiplied18(:,id_j) = dm_multiplied18(:,id_j) + v_ijk * dnormal_j + normal;
        dm_multiplied18(:,id_k) = dm_multiplied18(:,id_k) + v_ijk * dnormal_k + normal;
        % h (CoM multilied by mass)
        dh_multiplied24(:,id_i) = dh_multiplied24(:,id_i) + g .* dnormal_i + normal' .* dg(:,1:3);
        dh_multiplied24(:,id_j) = dh_multiplied24(:,id_j) + g .* dnormal_j + normal' .* dg(:,4:6);
        dh_multiplied24(:,id_k) = dh_multiplied24(:,id_k) + g .* dnormal_k + normal' .* dg(:,7:9);
        % rotational inertia
        % diagonal
        dDet_J_i = -vk*skew(vj);
        dDet_J_j = -vi*skew(vk);
        dDet_J_k = -vj*skew(vi);
        drot_g_diag_i = [dg(2,1:3) + dg(3,1:3); dg(1,1:3) + dg(3,1:3); dg(1,1:3) + dg(2,1:3)];
        drot_g_diag_j = [dg(2,4:6) + dg(3,4:6); dg(1,4:6) + dg(3,4:6); dg(1,4:6) + dg(2,4:6)];
        drot_g_diag_k = [dg(2,7:9) + dg(3,7:9); dg(1,7:9) + dg(3,7:9); dg(1,7:9) + dg(2,7:9)];
        drotInert_diag(:,id_i) = drotInert_diag(:,id_i) + rot_g_diag * dDet_J_i + det_J * drot_g_diag_i;
        drotInert_diag(:,id_j) = drotInert_diag(:,id_j) + rot_g_diag * dDet_J_j + det_J * drot_g_diag_j;
        drotInert_diag(:,id_k) = drotInert_diag(:,id_k) + rot_g_diag * dDet_J_k + det_J * drot_g_diag_k;
        % off-diagonal
        drotInert_off(:,id_i) = drotInert_off(:,id_i) - (r_ijk * dDet_J_i + det_J * r(v_ijk + vi));
        drotInert_off(:,id_j) = drotInert_off(:,id_j) - (r_ijk * dDet_J_j + det_J * r(v_ijk + vj));
        drotInert_off(:,id_k) = drotInert_off(:,id_k) - (r_ijk * dDet_J_k + det_J * r(v_ijk + vk));
    end
end
% phi = 10-dim vector
phi = [m_multiplied18/18; h_multiplied24/24; rotInert_diag/60; rotInert_off/120] * density;
if nargout >= 2
    dPhidVertex = [dm_multiplied18/18; dh_multiplied24/24; drotInert_diag/60; drotInert_off/120] * density;
end
end

function [g, dg] = g_vertex(vi, vj, vk)
if any(size(vi)~=[1 3]) || any(size(vj)~=[1 3]) || any(size(vk)~=[1 3])
    error(['g_vertex(): wrong input. size(v1)=' num2str(size(vi)) ', size(v2)=' num2str(size(vj)) ', size(v3)=' num2str(size(vk))])
end
g = transpose(vi.*vi + vi.*vj + vi.*vk + vj.*vj + vj.*vk + vk.*vk);
if nargout == 2
    % dg (g differentiated by vertices). size = 3 x 9
    v_ijk = vi+vj+vk;
    dg = [diag(v_ijk+vi) diag(v_ijk+vj) diag(v_ijk+vk)];
end
end

function rmat = r(v)
if any(size(v)~=[1 3])
    error('Error 1')
end
rmat = [v(2) v(1) 0;
    v(3) 0 v(1);
    0 v(3) v(2)];
end