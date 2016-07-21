function [canal_quadriga, canal_idd,ch_par] = generaCanalQuadriga(precision,ntx,nrx,realizaciones,random,output,los)

addpath('D:\TFM\Quadriga_new\source');  % Hay que añadir el path donde se encuentren los ficheros de Quadriga

s = simulation_parameters;
s.drifting_precision = precision; %Supondremos ondas esfericas en Rx y TX de forma
% que cada elemento del array posee unos angulos de salida y llegada
% diferentes

l = layout(s);

aTX = array; % Definición del array
aTX.generate('omni',ntx); % Generación de un array compuesto por antenas
% omnidirecionales de 30 elementos radiantes.

% Se deben de espacia equiespacial cada elemento (distancia entre elementos
% 0.5lambda

aTX.element_position(2,:) = [(-ntx+1)/2*s.wavelength:s.wavelength:ntx/2*s.wavelength];

%aTX.element_position(2,:) = linspace(-10,10,ntx);

% Se añade la antena TX al layout:
l.tx_array = aTX;
l.tx_position = [0; 0; 25];  % 25 metros de altura de la BS

% Ahora vamos con los receptores, se define un numero entero de MTs y se
% colocan de forma aleatoria en la celda:
l.simpar.show_progress_bars = 0;
l.no_rx = nrx;
% l.rx_array = aTX; l.rx_array.no_elements = 2;
% l.rx_array.element_position(2,:) = linspace(-5,5,2);
l.rx_position(3,:) = 0; % Colocamos todos los RX con altura 0;


if random
    l.randomize_rx_positions(2100,600,0,0,0);
else
    dx = randi(100,1);
    dy = randi(50,1);
    for y = 1:nrx
        
        l.rx_position(1,y) = dx+randi(10);
        l.rx_position(2,y) = dy+randi(10);
    end
end
l.rx_position(3,:) = 0; % Colocamos todos los RX con altura 0;
% A continuación se va a definir el escenario para poder, a partir de
% los parametros que se derivan de las mediciones, realizar la
% simulación de la comunicación.

for q = 1:size(l.track,2)
    if strcmp(los,'LOS')
        l.track(q).no_snapshots = 1;
        l.track(q).scenario = 'WINNER_UMi_B1_LOS';
    else
        l.track(q).scenario = 'WINNER_UMi_B1_NLOS';
        l.track(q).no_snapshots = 1;
        
    end
end

% Para una posición fija de los usuarios realizamos los calculos de la
% capacidad generando los coeifientes del canal en cada iteración
[a,ch_par] = l.generate_parameters;

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

for j=1:realizaciones
    
    if output
        %Output the simulation progress
        disp(['Realización = ' num2str(j)]);
    end
    % Obtenemos los LSPs y la instancia de la clase channel_builder:
    [h_channel, h_cb] = ch_par.get_channels;
    
    % Generamos tambien la mtriz de canal H para el caso iid Rayleigh para
    % comparar resultados con Quadriga
    H_idd_ideal =  (randn(ntx,nrx)+(sqrt(-1))*randn(ntx,nrx))/sqrt(2);  % parte de small scale fading
    
    % modelada como Rayleigh
    
    %% Se calcula la matriz de canal para esta iteración
    for k=1:l.no_rx
        h_user_paths = h_channel(1,k).coeff; % se seleciona el canal para cada usuario
        h_user = sum(h_user_paths,3); % se suman las contribuciones de los path
        H(:,k) = h_user;
        H_idd(:,k) = sqrt(l_scale(k)).*H_idd_ideal(:,k);
    end
    canal_quadriga{j} = H;
    canal_idd{j}=H_idd;
    % Las matrices de canal son de tamaño M x K;
end