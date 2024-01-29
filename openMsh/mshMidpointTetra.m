function [meshr,Ir] = mshMidpointTetra(mesh,I)

if nargin == 1
   I = (1:size(mesh, 1))'; 
end

% Check dimenion
if (size(mesh,2) ~= 4)
    error('mshMidpoint : unavailable case 1')
end

% Save color and replace by hierarchy
col      = mesh.col;
mesh.col = (1:length(mesh))';

% Submeshing 
meshs = mesh.sub(I);

% Edge meshes
edgs       = meshs.edg;
[edg,Itri] = mesh.edg;

% % Interface with edge multiplicity for triangle
[int,Iedg] = intersect(edg,edgs);
ind        = ismember(Itri,Iedg);
mlt        = sum(ind,2);

% 
% % Security
% tmp = setdiff(int,edgs);
% if (size(tmp.elt,1) ~= 0)
%     error('mshMidpoint.m : unavailable case 2');
% end

% Initialize refined mesh for element without refinement
meshr = mesh.sub(mlt==0);

% % Subdivision with 1 common edge
% tmp   = mesh.sub(mlt==1);
% tmp   = mshMidpoint1(tmp,int);
% meshr = union(meshr,tmp);
% 
% % Subdivision with 2 common edges
% tmp   = mesh.sub(mlt==2);
% tmp   = mshMidpoint2(tmp,int);
% meshr = union(meshr,tmp);
% 
% % Subdivision with 3 common edges
% tmp   = mesh.sub(mlt==3);
% tmp   = mshMidpoint3(tmp);
% meshr = union(meshr,tmp);

tmp = mesh;
tmp = mshMidpoint_tetra(tmp);
meshr = union(meshr,tmp);


% Parent indices and replace colours
Ir        = meshr.col;
meshr.col = col(Ir);

% Security
if (sum(mesh.ndv)-sum(meshr.ndv))/sum(mesh.ndv) > 1e-15*length(meshr)
    error('mshMidpoint.m : unavailable case 3');
end
end






function mesh = mshMidpoint_tetra(mesh)    
% Mesh nodes and edges center
[nds,ctr] = data(mesh);
   
% Refined mesh initialization with centered octahedron (4 tetrahedra)
Nelt = length(mesh);
elt  = reshape((1:4*Nelt)',Nelt,4);
col  = mesh.col;

% mesh_old = mesh;

vtx1 = [nds{1}; ctr{6}; ctr{4}; ctr{5}];
% tmp1 = msh(vtx1,elt,col);
mesh = msh(vtx1,elt,col);

vtx2 = [nds{2}; ctr{1}; ctr{3}; ctr{6}];
tmp2 = msh(vtx2,elt,col);
mesh = union(mesh,tmp2);

vtx3 = [nds{3}; ctr{2}; ctr{4}; ctr{1}];
tmp3 = msh(vtx3,elt,col);
mesh = union(mesh,tmp3);

vtx4 = [nds{4}; ctr{5}; ctr{3}; ctr{2}];
tmp4 = msh(vtx4,elt,col);
mesh = union(mesh,tmp4);

% vtx5 = [ctr{1}; ctr{2}; ctr{3}; ctr{5}];
% tmp5 = msh(vtx5,elt,col);
% mesh = union(mesh,tmp5);
% 
% vtx6 = [ctr{2}; ctr{5}; ctr{3}; ctr{6}];
% tmp6 = msh(vtx6,elt,col);
% mesh = union(mesh,tmp6);
% 
% vtx7 = [ctr{4}; ctr{6}; ctr{1}; ctr{5}];
% tmp7 = msh(vtx7,elt,col);
% mesh = union(mesh,tmp7);
% 
% vtx8 = [ctr{6}; ctr{1}; ctr{4}; ctr{2}];
% tmp8 = msh(vtx8,elt,col);
% mesh = union(mesh,tmp8);
vtx5 = [ctr{1}; ctr{3}; ctr{4}; ctr{2}];
tmp5 = msh(vtx5,elt,col);
mesh = union(mesh,tmp5);

vtx6 = [ctr{5}; ctr{4}; ctr{3}; ctr{2}];
tmp6 = msh(vtx6,elt,col);
mesh = union(mesh,tmp6);

vtx7 = [ctr{1}; ctr{3}; ctr{4}; ctr{6}];
tmp7 = msh(vtx7,elt,col);
mesh = union(mesh,tmp7);

vtx8 = [ctr{5}; ctr{4}; ctr{3}; ctr{6}];
tmp8 = msh(vtx8,elt,col);
mesh = union(mesh,tmp8);
% Security
if size(mesh.elt,1) ~= 8*Nelt
    error('mshMidpoint_tetra.m : unavailable case 1')
end
end






function [nds,ctr] = data(mesh)
% Triangle nodes
nds{1} = mesh.vtx(mesh.elt(:,1),:);
nds{2} = mesh.vtx(mesh.elt(:,2),:);
nds{3} = mesh.vtx(mesh.elt(:,3),:);
nds{4} = mesh.vtx(mesh.elt(:,4),:);

% Edges center
ctr{1} = 0.5 * (nds{2} + nds{3});
ctr{2} = 0.5 * (nds{3} + nds{4});
ctr{3} = 0.5 * (nds{4} + nds{2});
ctr{4} = 0.5 * (nds{3} + nds{1});
ctr{5} = 0.5 * (nds{1} + nds{4});
ctr{6} = 0.5 * (nds{1} + nds{2});

end

