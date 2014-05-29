%% Datos
clear;
close all;


% Se cargan las esperanzas, las varianzas y las covarianzas de las 3
% acciones
mu = [0.0427; 0.0015; 0.0285];
sigma = [0.1^2 0.0018 0.0011; 0.0018 0.1044^2 0.0026; 0.0011 0.0026 0.1411^2];

%% Por lagrange

Alag = [2*sigma ones(size(sigma,1),1); ones(1, size(sigma,2)) 0];
blag = [zeros(size(sigma,1),1); 1];
xlag = Alag\blag

%% Matrix Lagrange

auxmatlag = sigma\ones(size(sigma,2),1);

x0 = rand(size(auxmatlag));
%Resolver por householder, jacobi y gauss-seidel
%auxmatlag = householder(sigma, ones(size(sigma,2),1));
%auxmatlag = jacobi(sigma, ones(size(sigma,2),1), x0);
auxmatlag = gausssei(sigma, ones(size(sigma,2),1),x0);

xmatlag = auxmatlag/(ones(1,size(sigma,1))*auxmatlag)

%% Por parametrizacion

Z = [1 1; -1 0; 0 -1];
nuevaQ = Z'*sigma*Z;
bparam = -1*([1 zeros(1,size(sigma,1)-1)]*sigma*Z);

yparam = nuevaQ\bparam';
x0 = rand(size(nuevaQ),1);
%Resolver por householder, jacobi y gauss-seidel
%yparam = householder(nuevaQ,bparam');
%yparam = jacobi(nuevaQ,bparam',x0);
yparam = gausssei(nuevaQ,bparam',x0);

xparam = Z*yparam+[1; 0; 0]

