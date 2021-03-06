function Lch = CholFromLU(A)
    % Conseguir factorizacion LU de A
    [L,U,P] = lu(A);
    
    D = diag(diag(U));
    % Conseguir L de cholesky a partir de LU
    % Considerar la funcion diag de MATLAB
    
    Dsqr = sqrt(D);
    
    Lch = L*Dsqr;
    
    
    
    % Codigo para chequar que dio bien
    Ach = chol(A,'lower');
    
    for i = size(A,1)
        for j = size(A,2)
            if abs(Ach(i,j)-Lch(i,j)) > 0,001;
                error('Cholesky mal hecho');
            end
        end
    end
end
