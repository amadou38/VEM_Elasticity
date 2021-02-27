% Main function
% AUTEUR : Diallo Amadou, 28/09/2020

clear; clc; close all;
addpath('meshes/');
addpath('maillages/');
addpath('../../Mesh/mesh_files/');
mesh_file = 'mesh1.mat';
mesh = load(mesh_file);
% Test constantes
% E1 = 126; E2 = 11; G = 6.6; nu1 = 0.28; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2;  % constants
%E1 = 11.6; E2 = 0.9; G = 0.75; nu1 = 0.37; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 390;  % constants Epicea
E1 = 13.7; E2 = 2.24; G = 1.61; nu1 = 0.3; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 750;  % constants Hetre
%E1 = 11.5; E2 = 0.47; G = 0.5; nu1 = 0.005; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 392;  % constants Epicea
%E1 = 8.86; E2 = 0.54; G = 1.6; nu1 = 0.005; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 691;   % constants Sapin

C = [E1 E1*nu2 0; E2*nu1 E2 0; 0 0 G*c]/c;
%C = [3 1 0; 1 3 0; 0 0 1];
%C = [128.942 5.253 0; 5.253 13.309 0; 0 0 6.6];
%C = [2870 1280 0; 1280 3560 0; 0 0 1330];

%mesh.vertices(:,2) = 1.1*mesh.vertices(:,2); 
[U,Uex,ndofs] = elasto1(C,mesh);               % Calcul des solutions U et Uex
U1 = U(1:ndofs); U2 = U(ndofs+1:2*ndofs);
U1ex = Uex(1:ndofs); U2ex = Uex(ndofs+1:2*ndofs);
E_inf = max(abs(U1 - U1ex))/max(abs(U1ex))   % Error E_inf
E_L2 = norm(U1-U1ex)/norm(U1ex)              % Error E_L2
%[EV,V,ev,ndofs] = elasto2(C,mesh);            % Recherche de vp et de Vp
% FEV = sqrt(EV)/(2*pi);
% Fev = rho*sqrt(ev)/(2*pi);
% Affichage des resultats
figure
plot_sol(mesh, U1ex);
title('Solution - U1exacte ');
figure
plot_sol(mesh, U1);
title('Solution - U1app ');
% figure
% plot_sol(mesh, U2ex);
% title('Solution - U2exacte ');
% figure
% plot_sol(mesh, U1);
% title('Solution - U1app ');
% figure
% plot_sol(mesh, U2);
% title('Solution - U2app ');