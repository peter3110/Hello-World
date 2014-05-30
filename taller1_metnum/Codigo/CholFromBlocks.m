function Lch = CholFromBlocks(A)
    % Obtengo cholesky de una dimension menos
    L = chol(A(1:size(A,1)-1,1:size(A,2)-1), 'lower');
    
    % Calculo los valores que faltan y armo el resultado
    ultFila = (inv(L)*A(1:(size(A,1)-1), size(A,2)))';
    %ultFila = A(size(A,1),(1:(size(A,2)-1)))/L';
    ultElem = sqrt(A(size(A,1),size(A,2)) - (ultFila*ultFila'));
    Lch = [L zeros(size(A,1)-1,1); ultFila ultElem];
    
    
    % Codigo para chequar que dio bien
    Ach = chol(A,'lower');
    for i = size(A,1)
        for j = size(A,2)
            if abs(Ach(i,j)-Lch(i,j)) > 0,0001
                error('Cholesky mal hecho')
            end
        end
    end
end

