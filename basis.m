function p1 = basis(Xe,he)
% Function calculant les fonctions de bases interieures: p1
% 
% SYNOPSIS: p1 = basis(Xe,he);
% INPUT   : Xe   : le centroide               .he  : diametre
% OUTPUT  : p1   : les fonctions de la base P1
% AUTEUR : Diallo Amadou, 28/09/2020

p1 = {@(x,y) 1, @(x,y) 0, @(x,y) -(y-Xe(2))/he, @(x,y) (y-Xe(2))/he, @(x,y) (x-Xe(1))/he, @(x,y) 0; ...
      @(x,y) 0, @(x,y) 1, @(x,y)  (x-Xe(1))/he, @(x,y) (x-Xe(1))/he, @(x,y) 0, @(x,y) (y-Xe(2))/he };
  
end