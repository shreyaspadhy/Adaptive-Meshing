function [vert, idx, ne, nh, mua_node, mus_node, Q, mvec] = load_mesh_data(mesh, qpos, mpos, mua, mus, basis, refind)

nh = mesh.NodeCount();
ne = mesh.ElementCount();

[vert, idx] = mesh.Data();

mua_node = basis.Map('B->M', mua);
mus_node = basis.Map('B->M', mus);

mesh.SetQM(qpos, mpos);

Q = real(mesh.Qvec('Neumann','Gaussian',2));
mvec = real(mesh.Mvec('Gaussian',2,refind));

end
