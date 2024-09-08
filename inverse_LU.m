function A_inv = inverse_LU(A)
    [L, U, P] = lu(A); % Factorización LU con pivoteo parcial
    
    % Inversa de la matriz utilizando el operador \
    A_inv = zeros(size(A)); % Inicializar matriz inversa
    for i = 1:size(A, 1)
        e = zeros(size(A, 1), 1);
        e(i) = 1; % Vector de la i-ésima columna de la matriz identidad
        
        % Resolver sistema de ecuaciones lineales y asignar a la columna i de A_inv
        A_inv(:, i) = U \ (L \ (P * e));
    end
end