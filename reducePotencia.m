%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      Script: reducePotencia                       %%%%
%%%%                                                                   %%%%
%%%%        Script en el que se muestra el rendimiento de un sistema   %%%%
%%%%        MIMO masivo cuando se realiza el escalado de la potencia   %%%%
%%%%        tanto en downlink como en uplink. Se puede comprobaro      %%%%
%%%%        como el rendimiendo permenece constante aun que se este    %%%%
%%%%        realizando un escalado en función del número de antenas    %%%%
%%%%        en la estación base en uplink y escalado por M y           %%%%
%%%%        multiplicado por K en downlink.                            %%%%
%%%%                                                                   %%%%
%%%%           En este script se utiliza la normalización de           %%%%
%%%%           de las matrices de canal para que solo contengan        %%%%
%%%%               los efectos de small scale fading                   %%%%
%%%%  No se calculan los resultados para DPC (alto tiempo de ejecución)%%%%  
%%%%                                                                   %%%%
%%%%                  Resultados mostrados en 7.4.3                    %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
clc

addpath('D:\TFM\Quadriga')

% Parametros
nrx = 8;
precision = 5;
iteraciones = 2;
realizaciones = 10;
SNRdB =5;
M = [10:50:310]; %Vector de M
DL = 1;

for m = 1:length(M)
    
    for i = 1:iteraciones
        
        % Generamos los canales para una SNR dada (varias realizaciones donde
        % los usuarios se mantienen en la misma posición (large scale fading
        % iguales) pero small scale fading diferentes:
        
        [H_struct,H_iid_struct,x] = generaCanalQuadriga(precision,M(m),nrx,realizaciones,1,0,0);
        for h = 1 : realizaciones
            H_iid_struct{h} = (randn(M(m),nrx)+1i*randn(M(m),nrx))/sqrt(2);
        end
        % Con las realizaciones del canal se va a realizar la normalización
        % de todas ellas para conseguir que la potencia media de cada
        % usuario en todas ellas sea igual a 1.
        
        [canales_normalizados] = realizaNormalizacion(H_struct);
        [canales_normalizados_iid] = realizaNormalizacion(H_iid_struct);
        
        % Debido al numero tan alto realizaciones se van a coger 10 al azar
        % que seran las utilizadas para obtener la eficiencia espectral;
        
        index = randperm(realizaciones,10);
        % Calcular el sum rate para cada realización del canal y realizar la
        % media:
        
        for k = 1:length(index)
            
            H = canales_normalizados{index(k)};
            H_iid = canales_normalizados_iid{index(k)};

            if DL == 1
                SNRlin = 10^(SNRdB/10)*nrx/M(m); % Reducción de la potencia
                SNRdB_reducido = 10*log10(SNRlin);
                [C_mrt(k),C_zf(k),C_mmse(k)] = sum_rate_DL_sinDPC(H,SNRdB_reducido);
                [C_mrt_iid(k),C_zf_iid(k),C_mmse_iid(k)] = sum_rate_DL_sinDPC(H_iid,SNRdB_reducido);
            else
                SNRlin = 10^(SNRdB/10)/M(m); % Reducción de la potencia
                SNRdB_reducido = 10*log10(SNRlin);
                [C_mrt(k),C_zf(k),C_mmse(k),C_dpc(k)] = sum_rate_UL(H,SNRdB_reducido);
                [C_mrt_iid(k),C_zf_iid(k),C_mmse_iid(k),C_dpc_iid(k)] = sum_rate_UL(H_iid,SNRdB_reducido);
            end
        end
        
        CMRT(i) = mean(C_mrt);
        CZF(i) = mean(C_zf);
        CMMSE(i) = mean(C_mmse);
        
        CMRT_iid(i) = mean(C_mrt_iid);
        CZF_iid(i) = mean(C_zf_iid);
        CMMSE_iid(i) = mean(C_mmse_iid);
        
        
    end
    
    mrt(m) = mean(CMRT);
    zf(m) = mean(CZF);
    mmse(m) = mean(CMMSE);
    
    mrt_iid(m) = mean(CMRT_iid);
    zf_iid(m) = mean(CZF_iid);
    mmse_iid(m) = mean(CMMSE_iid);

end


figure
plot(M,mrt,'--o','LineWidth',1.5)
hold on
plot(M,zf,'o--r','LineWidth',1.5)
plot(M,mmse,'o--k','LineWidth',1.5);
grid on;
axis([0 max(M) 0 max(mmse)+5])
xlabel('M')
ylabel('Eficiencia espectral (bits/s/Hz)');
if DL ==1
    legend('MRT','ZF','MMSE');
    title('Eficiencia espectral con SNR variante en función de M (downlink)')
else
    legend('MRC','ZF','MMSE');
    title('Eficiencia espectral con SNR variante en función de M (uplink)')
end
