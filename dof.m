function D = dof(ne,np,p,Verts)
% Function calculant la matrice des degres de liberte: D_ra = dof_r(p_a)
% 
% SYNOPSIS: D = dof(ne,np,p,Verts);
% INPUT   : ne  : nbre de noeuds    .np   : nbre de monomes  
%           p  : base de monomes    .Verts: coordonees (x,y) de l element E              
%           Xe   : le centroide
% OUTPUT  : D   : Matrice de degre de liberte
% AUTEUR : Diallo Amadou, 28/09/2020

D = zeros(2*ne,2*np);

for k = 1:ne
    Vtx = Verts(k,:); % 
    for i = 1:2
        for j = 1:2*np
            D(2*k-2+i,j) = p{i,j}(Vtx(1),Vtx(2));
        end
    end
end

end

