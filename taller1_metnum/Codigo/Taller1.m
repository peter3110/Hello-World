%% leo imagen de entrada
close all, clear all;
% I = imread('../Senales/64x64/test.pgm');
% I = imread('../Senales/64x64/blond.pgm');
% I = imread('../Senales/64x64/canalet1.pgm');
% I = imread('../Senales/128x128/camera.pgm');
% I = imread('../Senales/256x256/lena.pgm');
% I = imread('../Senales/256x256/barbara.pgm');
figure,imshow(I),title('Imagen original');
%% agrego ruido uniforme 
sr = 64; % a partir de 64 (aprox) el filtrado vale la pena
IR = double(I) + randi([-sr,sr], size(I));
IR = uint8(IR);
figure,imshow(IR),title(['Imagen ruidosa. Interv=+-' num2str(sr)]);
%% Armado del sistema
Lambda = 1;
Dim = size(IR);
NInc = prod(Dim);


% Armado del    vector resultado
Utilde = Lambda*double(IR(:));

% Armado de la matriz a resolver
B = -1*ones(NInc,5);
B(:,3) = (Lambda+4)*ones(NInc,1);
d = [-Dim(1) -1 0 1 Dim(1)];
A = spdiags(B,d,NInc,NInc);


%% Vectores aleatorios para ver si es defpos (solo para imagenes de 64x64)

tic
for i = 0:100
    x = rand(size(A,2),1);
    if x ~= 0 
        if(x'*A*x <= 0)
            error('No Definida Positiva');
        end
    end        
end
toc       

%% Sylvester's Criterion (solo para imagenes de 64x64 y aun asi puede tardar unos minutos)

tic
for i = 1:size(A,1)
     if(det(A(1:i,1:i))<=0)
         error('error')
     end
end
toc



% %% Indices de perfil (solo para imagenes de 64x64)
% 
tic
L = chol(A,'lower');
for i = 1:size(L,1)
    if find(abs(L(i,:))>0.001,1) ~= find(abs(A(i,:))>0.001,1)
        error('Error en el indice de perfil')
    % A rellenar
    end
end
toc

%% Resolucion del sistema




% Resolucion directa
%

% Resolucion cholesky de matlab
%L = chol(A,'lower');
%Usol = (L*L')\Utilde;
 
% Resolucion cholesky por bloques

L = CholFromBlocks(A);
Usol = (L*L')\Utilde;


% Resolucion cholesky LU
%L = CholFromLU(A);
%Usol = (L*L')\Utilde;   





% Escalamiento y reshape para mostrar la imagen
Usol = Usol - min(Usol);
Usol = Usol / max(Usol);
Usol = uint8(Usol*255);
IFS = reshape(Usol,Dim);




figure,imshow(IFS);
figure,imshow([I,IR,IFS]),title(['Imagen filtrada, lambda=' num2str(Lambda)]);
%% calculo psnr
[p,m] = psnr(I(2:end-1,2:end-1), IR(2:end-1,2:end-1));

fprintf('Calidad imagen ruidosa.\n');
fprintf('PSNR=%g, ECM=%g\n',p,m);

[p,m] = psnr(I(2:end-1,2:end-1), IFS(2:end-1,2:end-1));
fprintf('Calidad imagen filtrada.\n');
fprintf('PSNR=%g, ECM=%g\n',p,m);


%%
figure;
S1 = edge(I,'canny');
S2 = edge(IR,'canny');
S3 = edge(IFS,'canny');
imshow([S1,S2,S3]);
