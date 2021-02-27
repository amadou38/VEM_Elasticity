function ML2 = ML2elt(np,p,ne,Verts,Xe,Area)
% Fonction qui calcule la matrice de masse elementaire
%
% SYNOPSIS: ML2 = ML2elt(np,p,ne,Verts,Xe,Area);
% INPUT   : ne   : nbre de noeuds        .Area: aire de E  
%           Verts: coordonees (x,y) de E .p: base des monomes            
%           Xe   : le centroide                     .np  : nbre de monomes
% OUTPUT  : ML2    : Matrice de masse
% AUTEUR : Diallo Amadou, 28/09/2020

RHS = RHS_P0(np,p,ne,Verts,Xe);
%RHS = RHS_P(C,ne,np,p,Verts,Xe,he);
%H = LHS(C,ne,np,p,Verts,Xe,he);
D = dof(ne,np,p,Verts);
H = LHS_P0(np,p,Verts,Xe,Area);
P0 = H\RHS;

ML2 = P0'*H*P0 + Area*(eye(2*ne) - D*P0)'*(eye(2*ne) - D*P0);

end