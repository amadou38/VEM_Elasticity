function B = RHS_P(C,ne,np,p,Verts,he)
% Fonction elementaire qui calcule le membre de droite de la matrice de
% projection P
% 
% SYNOPSIS: B = RHS_P(C,ne,np,p,Verts,he);
% INPUT   : ne  : nombre de points du polygone
%           np    : nombre de polynomes
%           Verts : polygone
% OUTPUT  : B   : la matrice
% AUTEUR : Diallo Amadou, 28/09/2020

B = zeros(2*np,2*ne);
Gs = zeros(3,2*np);
Gs(:,1) = [0 0 0]'; Gs(:,2) = Gs(:,1); Gs(:,3) = Gs(:,1); Gs(:,4) = 2*[0 0 1]'/he; 
Gs(:,5) = [1 0 0]'/he; Gs(:,6) = [0 1 0]'/he; 
Gs = C*Gs;
wrap = @(x,y) mod(x-1,y) + 1;
for j = 1:3
    for k = 1:ne
        B(j,2*k-1) = p{1,j}(Verts(k,1),Verts(k,2))/ne;
        B(j,2*k) = p{2,j}(Verts(k,1),Verts(k,2))/ne;
    end
end
for k = 1:ne
    
    Next = Verts(wrap(k+1, ne), :);
    Prev = Verts(wrap(k-1, ne), :);
    Avg = norm(Prev-Next); % Average
    Vec_normal = [Next(2) - Prev(2), Prev(1) - Next(1)]/Avg; % Normal vect
    Mat_normal = [Vec_normal(1) 0 Vec_normal(2); 0 Vec_normal(2) Vec_normal(1)]';
    for j = 4:2*np
        B(j, 2*k-1) = B(j, 2*k-1) + [0.5, 0]*(Avg*(Gs(:,j)'*Mat_normal)');
        B(j, 2*k) = B(j, 2*k) + [0, 0.5]*(Avg*(Gs(:,j)'*Mat_normal)');

    end
end

end