function M0 = Consis_term(C,ne,np,p,Verts,Xe,he)
% Function calculant le terme de consistance pour l elasticite orthotrope
% 2D: M0
% SYNOPSIS: M0 = Consis_term(C,ne,np,p,Verts,Xe,he);
% INPUT   : C    : la matrice des deformations        .ne  : nbre de noeuds
%           np   : nbre de monomes                    .p  : base de monomes
%           Verts: coordonees (x,y) de l element E              
%           Xe   : le centroide                       .he  : diametre
% OUTPUT  : M0   : Matrice de consistance
% AUTEUR : Diallo Amadou, 28/09/2020

Proj = projection(C,ne,np,p,Verts,Xe,he);
G = LHS(C,ne,np,p,Verts,Xe,he);
G(1:3,:) = 0;

M0 = Proj'*G*Proj;

end