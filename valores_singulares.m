%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Script: CalculaSingularSpread                      %%%% 
%%%%                                                                   %%%%
%%%%        Script en el que se calcula el numero de condicion         %%%%
%%%%        para diferentes configuraciones de MIMO multi-usuario      %%%%
%%%%        Para cada una de ellas, se realiza la normalización ya     %%%%
%%%%      que sólo interesa conocer información sobre la ortogonalidad %%%%
%%%%        de los usuarios y esto lo propocina la matriz H.           %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%            Resultados mostrados en sección: 7.4.2                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Parametros:
precision = 5;
v_ntx = [10 32 64 128];
nrx = 10;
realizaciones = 20;
random = 1;
output = 0;
los = 0;
iteraciones = 20;

for n = 1:length(v_ntx)
    
    ntx = v_ntx(n)
    
    for i = 1 : iteraciones
        
        % Generamos los canales para una SNR dada (varias realizaciones donde
        % los usuarios se mantienen en la misma posición (large scale fading
        % iguales) pero small scale fading diferentes:
        
        [H_struct,H_iid_struct,x,a] = generaCanalQuadriga(precision,ntx,nrx,realizaciones,1,0,0);
        
        for h = 1 : realizaciones
            H_iid_struct{h} = (randn(ntx,nrx)+1i*randn(ntx,nrx))/sqrt(2);
        end
        
        % Con las realizaciones del canal se va a realizar la normalización de
        % todas ellas para conseguir que la potencia media de cada usuario en
        % todas ellas sea igual a 1.
        
        [canales_normalizados] = realizaNormalizacion(H_struct);
        [canales_normalizados_iid] = realizaNormalizacion(H_iid_struct);
        
        % Debido al numero tan alto realizaciones se van a coger 5 al azar que
        % seran las utilizadas para obtener la eficiencia espectral;
        
        index = randperm(realizaciones,10);
        
        for m = 1:length(index)
            canales_quadriga{m} = canales_normalizados{index(m)};
            canales_iid{m} = canales_normalizados_iid{index(m)};
        end
        % Para cada realización se obtiene la distancia de Fraunhofer
        
        [condition_number(i), condition_number_idd(i)]= get_condition_number(canales_quadriga,canales_iid);
        
        eval(sprintf('condition_n_%d(i) =(condition_number(i));', ntx));
        eval(sprintf('condition_n_iid_%d(i) = (condition_number_idd(i));', ntx));
    end
    
    
end

% Medias y rango iqr:

mean_128 = median(condition_n_128);
mean_64 = median(condition_n_64);
mean_32 = median(condition_n_32);
mean_10 = median(condition_n_10);

iqr_128 = iqr(condition_n_128);
iqr_64 = iqr(condition_n_64);
iqr_32 = iqr(condition_n_32);
iqr_10 = iqr(condition_n_10);

mean_128_iid = mean(condition_n_iid_128);
mean_10_iid = mean(condition_n_iid_10);
iqr_128_iid = iqr(condition_n_iid_128);
iqr_10_iid = iqr(condition_n_iid_10);

% Calculo de las CDFs:

[cdf,intervalos_128] = get_CDF(condition_n_128,iteraciones);
[cdf,intervalos_128_iid] = get_CDF(condition_n_iid_128,iteraciones);

[cdf,intervalos_64] = get_CDF(condition_n_64,iteraciones);
[cdf,intervalos_64_iid] = get_CDF(condition_n_iid_64,iteraciones);

[cdf,intervalos_32] = get_CDF(condition_n_32,iteraciones);
[cdf,intervalos_32_iid] = get_CDF(condition_n_iid_32,iteraciones);

[cdf,intervalos_10] = get_CDF(condition_n_10,iteraciones);
[cdf,intervalos_10_iid] = get_CDF(condition_n_iid_10,iteraciones);


figure
plot(intervalos_128,cdf,'Color',[1,0.7,0.1],'LineWidth', 1.5)
hold on
plot(intervalos_64,cdf,'--r','LineWidth', 1.5)
plot(intervalos_32,cdf,'x-k','LineWidth', 1.5)
plot(intervalos_10,cdf,'LineWidth', 1.5);
plot(intervalos_128_iid,cdf,'--','Color',[1,0.7,0.1],'LineWidth', 1.2)
plot(intervalos_64_iid,cdf,'--r','LineWidth', 1.2)
plot(intervalos_32_iid,cdf,'k--','LineWidth', 1.2)
plot(intervalos_10_iid,cdf,'--','LineWidth', 1.2)
axis([0 max(intervalos_10+2) 0 1])
grid on
title('CDF del número de condición en un escenario NLOS')
legend('128 antenas','64 antenas','32 antenas','10 antenas')
ylabel('Prob(condition number) < y')
xlabel('Condition number (dB)')
