function G = LHS(C,ne,np,p,Verts,Xe,he)
% Function calculant le premier terme de la projection: G = B*D
% 
% SYNOPSIS: G = LHS(C,ne,np,p,Verts,Xe,he);
% INPUT   : C    : la matrice des deformations        .ne  : nbre de noeuds
%           np   : nbre de monomes                    .p  : base de monomes
%           Verts: coordonees (x,y) de l element E              
%           Xe   : le centroide                       .he  : diametre
% OUTPUT  : G    : la matrice
% AUTEUR : Diallo Amadou, 28/09/2020

D = dof(ne,np,p,Verts);
B = RHS_P(C,ne,np,p,Verts,he);
G = B*D;
V = zeros(2*ne,2*np);
for k = 1:ne
    Vtx = Verts(k,:); % 
    for i = 1:2
        for j = 1:2*np
            V(2*k-2+i,j) = p{i,j}(Vtx(1),Vtx(2));
        end
    end
end

for i = 1:3
    G(i,:) = sum(V(:,i).*V)/ne;
end

end