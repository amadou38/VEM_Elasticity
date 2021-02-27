function G = LHS_P0(np,p,Verts,Xe,Area)
% Function calculant le premier terme de la projection P0: 
% 
% SYNOPSIS: G = LHS_P0(ne,np,p,Verts,Xe,he);
% INPUT   : ne   : nbre de noeuds
%           np   : nbre de monomes                    .p  : base de monomes
%           Verts: coordonees (x,y) de l element E              
%           Xe   : le centroide                       .Area  : |E|
% OUTPUT  : G    : la matrice
% AUTEUR : Diallo Amadou, 28/09/2020

G = zeros(2*np);
for i = 1:2*np
    for j = 1:2*np
        f = @(x,y) p{1,i}(x,y).*p{1,j}(x,y) + p{2,i}(x,y).*p{2,j}(x,y);
        G(i,j) = integral_fct(f, Verts, Xe, 3)/Area;
    end
end
% D = dof(np,p,ne,Verts,Xe,Area);
% RHS = RHS_P(nu,lambda,np,p,ne,Verts,he);
% G = RHS*D;

end