function M1 = Stab_term(C,ne,np,p,Verts,Xe,he)
% Function calculant le terme de stabilite pour l elasticite orthotrope
% 2D: M1
% SYNOPSIS: M1 = Stab_term(C,ne,np,p,Verts,Xe,he);
% INPUT   : C    : la matrice des deformations        .ne  : nbre de noeuds
%           np   : nbre de monomes                    .p  : base de monomes
%           Verts: coordonees (x,y) de l element E              
%           Xe   : le centroide                       .he  : diametre
% OUTPUT  : M1   : Matrice de stabilite
% AUTEUR : Diallo Amadou, 28/09/2020

Proj = projection(C,ne,np,p,Verts,Xe,he);
D = dof(ne,np,p,Verts);

M1 = (eye(2*ne) - D*Proj)'*(eye(2*ne) - D*Proj);

end