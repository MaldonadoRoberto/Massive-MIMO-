%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Función: generaCanalQuadriga                      %%%%
%%%%                                                                   %%%%
%%%%      Función en la que se genera la matriz de Quadriga según los  %%%%
%%%%      parametros introducidos por el usuario. Además de obtener    %%%%
%%%%      el canal a traves de la simulación de Quadriga genera un     %%%%
%%%%      canal i.i.d Rayleigh al que se le incorporan los parametros  %%%%
%%%%      de large scale fading según explicó Marzetta en su paper:    %%%%
%%%%           ' Noncooperative Cellular Wireless with Unlimited       %%%%
%%%%              Numbers of Base Station Antennas'                    %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE ENTRADA:                                       %%%%
%%%%                                                                   %%%%
%%%%      precision: tipo de suposición sobre la pragación de ondas    %%%%
%%%%          5 = ondas esfericas en ambos extremos.                   %%%%
%%%%                                                                   %%%%
%%%%      ntx: número de antenas en la estación base                   %%%%
%%%%                                                                   %%%%
%%%%      nrx: número de usuarios en el sistema.                       %%%%
%%%%                                                                   %%%%
%%%%      realizaciones: número de realizaciones para una distribución %%%%
%%%%      de usuarios dada.                                            %%%%
%%%%                                                                   %%%%
%%%%      random: variable booleana que decide si los usuarios se      %%%%
%%%%      situan aleatoriamente sobre la celda o no.                   %%%%
%%%%                                                                   %%%%
%%%%      output: variable booleana que decide si se muestra por       %%%%
%%%%      pantalla la salida de la función                             %%%%
%%%%                                                                   %%%%
%%%%      los: variable booleana que decide si se simula un entorno    %%%%
%%%%            LOS o NLOS.                                            %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE SALIDA:                                        %%%%
%%%%                                                                   %%%%
%%%%      canal_quadriga: struct que contiene en cada posición una     %%%%
%%%%      matriz de tamaño [M x K] obtenida mediante Quadriga          %%%%
%%%%                                                                   %%%%
%%%%      canal_iid: struct que contiene en cada posición una          %%%%
%%%%      matriz de tamaño [M x K] obtenida mediante iid + large scale %%%%
%%%%                                                                   %%%%
%%%%      ch_par: parametros del canal obtenidos mediante Quadriga     %%%%
%%%%                                                                   %%%%
%%%%      h_channel: parametros del canal generado por Quadriga        %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [canal_quadriga, canal_idd,ch_par,h_channel] = generaCanalQuadriga(precision,ntx,nrx,realizaciones,random,output,los)

%% Definición de los parametros de la simulación

addpath('..') % Hay que añadir el path donde se encuentren los ficheros de Quadriga

s = simulation_parameters;
s.drifting_precision = precision;

l = layout(s);

% Definición del array en BS
aTX = array;
aTX.generate('omni',ntx);
% Equiespaciado entre elementos de lmabda/2
aTX.element_position(2,:) = [(-ntx+1)/2*s.wavelength:s.wavelength:ntx/2*s.wavelength];


% Se añade la antena TX al layout:
l.tx_array = aTX;
l.tx_position = [0; 0; 25];  % 25 metros de altura de la BS

% Añade una matriz de acoplamiento en el array de la BS.

% l.tx_array.coupling = generaMatrizCoupling(ntx);

% Definición de los usuarios
l.simpar.show_progress_bars = 0;
l.no_rx = nrx;
l.rx_position(3,:) = 0; %Colocamos todos los RX con altura 0;

if random
    l.randomize_rx_positions(1000,250,0,0,0);  % Generación de la posición de los usuarios entre
    % disntacias de 250 y 1000 metros.
else
    % Generación de los usuarios situados próximos entre si. El valor
    % espaciado puede variar según como se quiera que esten los usuarios
    % de proximos
    dx = 500;
    dy = 500;
    espaciado = 25;
    for y = 1:nrx
        l.rx_position(1,y) = dx+randi(espaciado);
        l.rx_position(2,y) = dy+randi(espaciado);
    end
end
l.rx_position(3,:) = 0; % Colocamos todos los RX con altura 0;

% A continuación se va a definir el escenario para poder, a partir de
% los parametros que se derivan de las mediciones, realizar la
% simulación de la comunicación.
for q = 1:size(l.track,2)
    if los == 1
        l.track(q).no_snapshots = 1;
        l.track(q).scenario = 'WINNER_UMa_C2_LOS';
    else
        l.track(q).scenario = 'WINNER_UMa_C2_NLOS';
        l.track(q).no_snapshots = 1;
        
    end
end

[a,ch_par] = l.generate_parameters;

%% Generación de los parametros de large scale según Marzetta
% Large scale fading para el caso del idd;
% Large scale fading segun Marzeta
path_loss_dB = ch_par.get_pl;  % Se obtienen los vlaores de path
% loss en dB

% Obtenemos el factor shadow fading para cada receptor:
sf = ch_par.sf;

% Pasamos el pathloss a unidades lineales
path_loss = 10.^(path_loss_dB/10);

% Calculo de large scale:
l_scale = sf ./ path_loss';

%% Obtención de los canales
rng('shuffle');
for j=1:realizaciones
    
    if output
        %Output the simulation progress
        disp(['Realización = ' num2str(j)]);
    end
    % Obtenemos los LSPs y la instancia de la clase channel_builder:
    [h_channel, h_cb] = ch_par.get_channels;
    
    % Generamos tambien la mtriz de canal H para el caso iid Rayleigh para
    % comparar resultados con Quadriga
    H_idd_ideal =  (randn(ntx,nrx)+1i*randn(ntx,nrx))/sqrt(2);  % parte de small scale fading
    D = diag(l_scale);
    
    %% Se calcula la matriz de canal para esta iteración
    for k=1:l.no_rx
        h_user_paths = h_channel(1,k).coeff; % se seleciona el canal para cada usuario
        h_user = sum(h_user_paths,3); % se suman las contribuciones de los path
        H(:,k) = h_user;
    end
    canal_quadriga{j} = H;
    canal_idd{j}=H_idd_ideal*D^(1/2);   % Generación de la matriz de canal según Marzetta
end