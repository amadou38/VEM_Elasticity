function [EV,V,ev,ndofs] = elasto2(C,mesh,ind)
% Function calculant les vp et Vp de l elasticite orth: EV, V
% 
% SYNOPSIS: [EV,V,ev,ndofs] = elasto2(C,mesh);
% INPUT   : C  : matrice des deformations    .mesh   : maillage du domaine
%           ind: nbre d elements a calcuer (vp et Vp)
%
% OUTPUT  : EV    : valeurs propres
%           V     : Vecteurs propres
%           ndofs: dimension
% AUTEUR : Diallo Amadou, 28/09/2020

ndofs = length(mesh.vertices);
np = 3; % dim(P1(E))
nd = 2*ndofs;           % pour k = 1, car nd = dim(V_h(E)) = 2*nvtx + 2*nvtx*(k-1) + k(k-1)
M0 = zeros(nd);         % matrix for the consistency term
M1 = zeros(nd);
ML2 = zeros(nd);
for l = 1:length(mesh.elements)
    vtx_id = mesh.elements{l};
    [Verts,Xe,ne,he,Area] = Polygon(mesh,l);
    p = basis(Xe,he);
    M0_E = Consis_term(C,ne,np,p,Verts,Xe,he);
    alphaE = 0.5*trace(M0_E);
    %alphaE = 4.25;
    M1_E = alphaE*Stab_term(C,ne,np,p,Verts,Xe,he);
    ML2_E = ML2elt(np,p,ne,Verts,Xe,Area);
    for i = 1:ne
        I = vtx_id(i); Ip = ndofs+I;
        for j = 1:ne
            J = vtx_id(j); Jp = J + ndofs;
            M0(I,J) = M0(I,J) + M0_E(2*i-1,2*j-1);
            M0(Ip,J) = M0(Ip,J) + M0_E(2*i,2*j-1);
            M0(I,Jp) = M0(I,Jp) + M0_E(2*i-1,2*j);
            M0(Ip,Jp) = M0(Ip,Jp) + M0_E(2*i,2*j);
            M1(I,J) = M1(I,J) + M1_E(2*i-1,2*j-1);
            M1(Ip,J) = M1(Ip,J) + M1_E(2*i,2*j-1);
            M1(I,Jp) = M1(I,Jp) + M1_E(2*i-1,2*j);
            M1(Ip,Jp) = M1(Ip,Jp) + M1_E(2*i,2*j);
            ML2(I,J) = ML2(I,J) + ML2_E(2*i-1,2*j-1);
            ML2(Ip,J) = ML2(Ip,J) + ML2_E(2*i,2*j-1);
            ML2(I,Jp) = ML2(I,Jp) + ML2_E(2*i-1,2*j);
            ML2(Ip,Jp) = ML2(Ip,Jp) + ML2_E(2*i,2*j);
        end
    end
end
M = M0 + M1;
idofs = ~ismember(1:ndofs, mesh.boundary); % noeuds qui ne sont pas sur le bord
idofs = find(idofs);
idofs2d =[idofs,ndofs+idofs]; 

M = M(idofs2d,idofs2d); 
ML2 = ML2(idofs2d,idofs2d);
spy(ML2);
figure, spy(M);
Neig = 15;
tic
[V,ev] = eig(M,ML2);
toc
[ev, indices] = sort(diag(ev));
V = V(:,indices);
%ev = 1./sqrt(ev/rho);
EV = zeros(Neig,1);
W = zeros(size(M,1),Neig);
k = 0;
for i = 1:size(M,1)
    if ((ev(i) >=1e-2) && (k<Neig))
        k = k + 1;
        EV(k) = ev(i);
        W(:,k) = V(:,i);
    elseif (k == Neig)
        break;
    end
end
EV = EV';



end