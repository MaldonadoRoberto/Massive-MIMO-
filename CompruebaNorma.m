%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      Script: CompruebaNorma                       %%%%
%%%%                                                                   %%%%
%%%%        Script en el que se comprueba que la matriz de canal       %%%%
%%%%        proporcionada por Quadriga se puede descomponer en las     %%%%
%%%%        componentes de large scale fading + small scale fading     %%%%   
%%%%        para ello, se obtiene la potencia promedio sobre N         %%%%
%%%%        realizaciones del canal. Despues se obtiene la matriz      %%%%
%%%%       que solo contiene las contribuciones del large scale fading %%%%
%%%%                  aplicando H = G*D^(-1/2)                         %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%                  Resultados mostrados en 7.2.1                    %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
clc

% Iniciación de los parametros
M = 100;
K = 5;
precision = 5;
realizaciones = 100;


% Realizaciones de la misma matriz de quadriga;

[H_q,H_i,ch_par,x] = generaCanalQuadriga(precision,M,K,realizaciones,1,0,0);

% Calculo de la potencia para cada realizacion de matriz;

for i = 1:realizaciones
    realizacion = H_q{i};  
    % obtrención de la matriz que contiene la potencia de cad ausuario; 
    pow(:,:,i) = abs(realizacion).^2;
end

pow_media = mean(pow,3);

% comparación con el path gain;
path_loss_dB = ch_par.get_pl; % path loss
sf = ch_par.sf;  % shadow fading

% Pasamos el pathloss a unidades lineales
path_loss = 10.^(path_loss_dB/10);

% Calculo de large scale:
l_scale = sf ./ path_loss';


% compracion de resultaods
fprintf('Comparación de los resultados para comprobar si la matriz de Quadriga contiene las contribuciones de large scale fading\n\n')
fprintf('Matriz de Quadriga: \n')
pow_media_media = mean(pow_media)
fprintf('Large Scale fading (obtenido segun Marzetta): \n')
l_scale


% la matriz D por tanto es la matriz didagonal de pow_media_media;
D = diag(pow_media_media);

% esta D sera utilizada para normalizar el canal que propocionar quadrigra
% (cada realizacion);

for k = 1:realizaciones
    H_sin_pg(:,:,k) = H_q{k}*D^(-1/2);
end

% vuelvo a realizar la comprobación de la potencia de esta matriz (sin
% pathgain)

for i = 1:realizaciones
    
    realizacion = H_sin_pg(:,:,i);
    
    % obtrención de la matriz que contiene la potencia de cad ausuario;
    
    pow_sinpg(:,:,i) = abs(realizacion).^2;
    aa(:,:,i) = channel_norm(H_q{i},1);
end

pow_media_sinpg = mean(pow_sinpg,3);

fprintf('Potencia media de cada posición de la matriz de canal sin path gain:\n')
media_total = mean(pow_media_sinpg)

