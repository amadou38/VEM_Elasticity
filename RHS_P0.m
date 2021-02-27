function RHS = RHS_P0(np,p,ne,Verts,Xe)
% Fonction elementaire qui calcule le membre de droite de la matrice de
% projection P0
% 
% SYNOPSIS: RHS = RHS_P0(np,ne,Verts,he);
% INPUT   : ne  : nombre de points du polygone
%           np    : nombre de polynomes
%           Verts : polygone
% OUTPUT  : RHS   : la matrice
% AUTEUR : Diallo Amadou, 28/09/2020

RHS = zeros(2*np,2*ne);

for k = 1:ne
    for i = 1:2*np
        RHS(i,2*k-1) = integral_fct(p{1,i}, Verts, Xe, 3);
        RHS(i,2*k) = integral_fct(p{2,i}, Verts, Xe, 3);
    end
end
end