function Proj = projection(C,ne,np,p,Verts,Xe,he)
% Fonction qui calcule la matrice de projection
%
% SYNOPSIS: Proj = projection(C,ne,np,p,Verts,Xe,he);
% INPUT   : C   : matrice des deformations     .ne: nbre de noeuds  
%           Verts: coordonees (x,y) de E .p: base des monomes            
%           Xe   : le centroide                     .np  : nbre de monomes
%           he : diametre
% OUTPUT  : Proj    : Matrice de projection
% AUTEUR : Diallo Amadou, 28/09/2020

% D = dof(ne,np,p,Verts);
B = RHS_P(C,ne,np,p,Verts,he);
G = LHS(C,ne,np,p,Verts,Xe,he);
Proj = G\B;

end