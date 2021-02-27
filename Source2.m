function F = Source2(C,ne,Xe,Area)
% AUTEUR : Diallo Amadou, 28/09/2020

% Uex(x,y) = (u1,u1) = (sin(pi*x)sin(pi*y), sin(pi*x)sin(pi*y))
F = zeros(2*ne,1);
x = Xe(1); y = Xe(2);
u1 = pi*pi*sin(pi*x).*sin(pi*y); u = pi*pi*cos(pi*x).*cos(pi*y);
f1 = u1*(C(1,1)+C(3,3)) - u*(C(1,2)+C(3,3));% + 2*sin(pi*x).*sin(pi*y);
f2 = u1*(C(2,2)+C(3,3)) - u*(C(1,2)+C(3,3));% + 2*sin(pi*x).*sin(pi*y);
for i = 1:ne
    F(2*i-1) = f1*Area/ne;
    F(2*i) = f2*Area/ne;
end

end