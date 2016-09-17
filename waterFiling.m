%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Funcion: waterFilling                            %%%%
%%%%                                                                   %%%%
%%%%       Funcion en la que se las potencias de un sistema MIMO,      %%%%
%%%%       puede usarse para un sistema MIMO single user o para realizar %%
%%%%       la maximización en water filling iterativo donde se aplica  %%%%
%%%%       los canales efectivos de cada usuario.                      %%%%
%%%%                                                                   %%%%
%%%%     El metodo utilizado esta descrito en el apartado 2.3.5        %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE ENTRADA:                                       %%%%
%%%%                                                                   %%%%
%%%%      H: matriz de canal                                           %%%%
%%%%                                                                   %%%%
%%%%      SNR : SNR en escala lineal                                   %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE SALIDA:                                        %%%%
%%%%                                                                   %%%%
%%%%      potencias_optimas: vector de potencias que                   %%%%
%%%%      maximizan la capacidad                                       %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function potencias_optimas = waterFiling(H,SNR)
r=rank(H);
p=1;
[U,D,V] = svd(H);
lambda = diag(D).^2;

while p<r
    
    inv_lambda = 1./lambda(1:r-p+1);
    mu = (1/(r-p+1)) * (1+(1/SNR)*sum(inv_lambda)); % Calculo del nivel de agua
    potencias = repmat(mu,(r-p+1),1) - (1./(SNR.*lambda(1:r-p+1))); % Asignación de potencias
    
    if (potencias(r-p+1) < 0)  % Si la ultima potencia asignada es menor que 0 se descarta y se vuelve a iterar.
        potencias(r-p+1) = 0;
        p = p+1;
    else
        p = r;
    end
    
end

% Se da formato a los resultados ( se añaden 0 en las potencias que no han
% sido asignadas)
potencias = potencias.';
potencias_optimas = zeros(1,r);
potencias_optimas(1:size(potencias,2)) = potencias;

end