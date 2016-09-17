%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Script: correlation_coeff                          %%%% 
%%%%                                                                   %%%%
%%%%        Script que obtiene el coeficiente de correlación           %%%%
%%%%       entre dos pares de usuarios. Para obtener el coeficiente de %%%%
%%%%       correlación medio, se realiza una media sobre diferentes    %%%%
%%%%              coeficientes de correlación calculados               %%%%
%%%%            para pares de usuario elegidos aleatoriamente.         %%%%
%%%%                                                                   %%%%
%%%%           El método utilizado se basa en el paper: Channel        %%%% 
%%%%               Measurements for Large Antenna Arrays.              %%%% 
%%%%        Esta función realiza el calculo en función del             %%%%
%%%%             número de antenas en la estación base.                %%%%
%%%%    Se define el coeficiente de correlacion entre vectores como:   %%%% 
%%%%                                                                   %%%%
%%%%            delta = norma(h_k^H*h_j)/(norma(h_k)*norma(h_j));      %%%%
%%%%                                                                   %%%%
%%%%            Resultados mostrados en sección: 7.4.2                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Iniciación de los parametros
n_BS = 10:30:220; % Vector de antenas en la estación base
nrx = 10; % Número de usuarios en la celda
repeticiones = 20; % Repeticiones que se van a realizar para la misma configuración

% En este primer bucle nn toma dos valores; si nn = 1 => Se simula un
% entorno LOS; si por el contrado nn = 2 => Se simula un entorno NLOS.
for nn = 1:2
    
    for i = 1:length(n_BS)  % Se recorre el vector de M
        
        % Se toman vectores aleatorios y se calculan sus coeficientes de
        % correlacion
        
        for j = 1:repeticiones
            
            if nn == 1
                [canalQ,canalIDD]=generaCanalQuadriga(5,n_BS(i),nrx,1,0,0,1);
            else
                [canalQ,canalIDD]=generaCanalQuadriga(5,n_BS(i),nrx,1,0,0,0);
            end
            
            H = canalQ{1}; % Seleción del canal
            H_iid = (randn(n_BS(i),nrx)+1i*randn(n_BS(i),nrx))./sqrt(2); % Generación de i.i.d
            
            % Generación aleatoria de dos números para calcular la
            % correlación entre ellos:
            for m = 1:20
                canal_a_estudiar = zeros(1,2);
                while canal_a_estudiar(1) == canal_a_estudiar(2)
                    canal_a_estudiar = randi(nrx,[1 2]);
                end
                
                % Obtención del coeficiente de correlación
                delta(m) = norm(H(:,canal_a_estudiar(1))'*H(:,canal_a_estudiar(2)))/(norm(H(:,canal_a_estudiar(1)))*norm(H(:,canal_a_estudiar(2))));
                delta_iid(m) = norm(H_iid(:,canal_a_estudiar(1))'*H_iid(:,canal_a_estudiar(2)))/(norm(H_iid(:,canal_a_estudiar(1)))*norm(H_iid(:,canal_a_estudiar(2))));
                
            end
            
            % Media sobre todas las parejas de canales de usuario
            correlacion_rep(j) = mean(delta);
            correlacion_rep_iid(j) = mean(delta_iid);
            
            
        end
        
        % Media para cada repeticion
        correlation_coeff(i) = mean(correlacion_rep);
        correlation_coeff_iid(i) = mean(correlacion_rep_iid);
    end
    
    % Valor del coeficiente de correlación para un numero de antenas M
    % dado.
    correlation(:,nn) = correlation_coeff;
end


% OPCIONAL: Interpolación sobre los datos para obtener una representación
% más lineal de los resultados.
n_bs_nuevo = 2:0.5:220;
correlation_inter = interp1(n_BS,correlation,n_bs_nuevo,'spline');
correlation_iid_inter = interp1(n_BS,correlation_coeff_iid,n_bs_nuevo,'spline');


figure
plot(n_bs_nuevo,correlation_inter(:,1),'--','LineWidth',1.5)
hold on
plot(n_bs_nuevo,correlation_inter(:,2),'--r','LineWidth',1.5)
plot(n_bs_nuevo,correlation_iid_inter,'k','LineWidth',1.5);
plot(n_BS,correlation(:,1),'o')
plot(n_BS,correlation(:,2),'r*')
plot(n_BS,correlation_coeff_iid,'ok')
grid on
legend('LOS','NLOS','IID')
xlabel('Número de antenas en transmisión')
ylabel('Coeficiente de correlación')
title('Comparativa del coficiente de correlación frente a M')
axis([0 max(n_bs_nuevo)+2 0 1])

figure
plot(n_BS,correlation(:,1),'--o','LineWidth',1.5)
hold on
plot(n_BS,correlation(:,2),'-*r','LineWidth',1.5)
plot(n_BS,correlation_coeff_iid,'o--k','LineWidth',1.5)
grid on
legend('LOS','NLOS','IID')
xlabel('Número de antenas en transmisión')
ylabel('Coeficiente de correlación')
title('Comparativa del coficiente de correlación frente a M (sin interpolación)')
axis([0 max(n_BS)+2 0 1])