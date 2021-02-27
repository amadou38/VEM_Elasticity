function [Verts,Xe,Nvtx,he,Area] = Polygon(mesh,l)
% AUTEUR : Diallo Amadou, 28/09/2020

i_vert = mesh.elements{l};
Verts = mesh.vertices(i_vert,:);
Nvtx = length(i_vert);
Area = 0.5*abs(sum(Verts(:,1).*Verts([2:end,1],2) - Verts([2:end,1],1).*Verts(:,2)));
Xe = sum((Verts + Verts([2:end,1],:)) .* repmat(Verts(:,1).*Verts([2:end,1],2) - Verts([2:end,1],1).*Verts(:,2),1,2));
Xe = Xe/(6*Area);
he = 0;
for i = 1:Nvtx-1
    for j = (i+1):Nvtx
        he = max(he, norm(Verts(i, :) - Verts(j, :)));
    end
end


end