
clear; close all; clc;
L2 = [0.00792603591790700;0.00147804712576971;0.000815316821111028;0.000166781640021099;8.48431318628270e-05];
h = [ 100 500 1000 5000 10000];
loglog(h, L2, 'r--o',h, 1./h, 'k-')
legend('Erreur L2','Ordre 1')
xlabel('NELEM'); ylabel('ErreurL2');
title('Norme L2')