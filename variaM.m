%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      Script: variaM                               %%%%
%%%%                                                                   %%%%
%%%%        Script el que se obtiene la eficiencia espectral           %%%%
%%%%       para un sistema de MIMO masivo en función del numero de     %%%%
%%%%     antenas en BS. Se puede utilizar tanto para el enlace uplink  %%%%
%%%%        como para downlink. La SNR y el número de usuarios         %%%%
%%%%        permanece constante a lo largo de la simulación.           %%%%
%%%%                                                                   %%%%
%%%%           En este script se utiliza la normalización de           %%%%
%%%%           de las matrices de canal para que solo contengan        %%%%
%%%%               los efectos de small scale fading                   %%%%
%%%%               
%%%%   El tiempo de ejecución cuando se usa DPC (downlink) es elevado. %%%%
%%%%                                                                   %%%%
%%%%                  Resultados mostrados en 7.4.1                    %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

addpath('..')
rng('shuffle'); %Initiate the random number generators with a random seed

% Graficas de la tesis:  (downlink)
ntx = [10:20:120];
nrx = 8;
precision = 5;
iteraciones = 2;
realizaciones = 50;
SNRdB = [5];
DL = 1;

tic
for m = 1:length(ntx)
        fprintf('iteracion %d de %d\n',m,length(ntx));
    for i = 1:iteraciones
        
        % Generamos los canales para una SNR dada (varias realizaciones donde
        % los usuarios se mantienen en la misma posición (large scale fading
        % iguales) pero small scale fading diferentes:
        
        [H_struct,H_iid_struct,x,a] = generaCanalQuadriga(precision,ntx(m),nrx,realizaciones,1,0,0);
        for h = 1 : realizaciones
            H_iid_struct{h} = (randn(ntx(m),nrx)+1i*randn(ntx(m),nrx))/sqrt(2);
        end
        % Con las realizaciones del canal se va a realizar la normalización
        % de todas ellas para conseguir que la potencia media de cada
        % usuario en todas ellas sea igual a 1.
        
        [canales_normalizados] = realizaNormalizacion(H_struct);
        [canales_normalizados_iid] = realizaNormalizacion(H_iid_struct);
        
        % Calcular el sum rate para cada realización del canal y realizar la
        % media:
        
        % Debido al numero tan alto realizaciones se van a coger 5 al azar
        % que seran las utilizadas para obtener la eficiencia espectral;
        
        index = randperm(realizaciones,10);
        for k = 1:length(index)
            
            H = canales_normalizados{index(k)};
            H_iid = canales_normalizados_iid{index(k)};
            
            if DL == 1
                [C_mrt(k),C_zf(k),C_mmse(k),C_dpc(k)] = sum_rate_DL(H,SNRdB);
                [C_mrt_iid(k),C_zf_iid(k),C_mmse_iid(k),C_dpc_iid(k)] = sum_rate_DL(H_iid,SNRdB);
            else
                [C_mrt(k),C_zf(k),C_mmse(k),C_dpc(k)] = sum_rate_UL(H,SNRdB);
                [C_mrt_iid(k),C_zf_iid(k),C_mmse_iid(k),C_dpc_iid(k)] = sum_rate_UL(H_iid,SNRdB);
            end
        end
        
        CMRT(i) = mean(C_mrt);
        CZF(i) = mean(C_zf);
        CMMSE(i) = mean(C_mmse);
        CDPC(i) = mean(C_dpc);
        
        CMRT_iid(i) = mean(C_mrt_iid);
        CZF_iid(i) = mean(C_zf_iid);
        CMMSE_iid(i) = mean(C_mmse_iid);
        CDPC_iid(i) = mean(C_dpc_iid);
        
    end
    
    CMRT_media(m) = mean(CMRT);
    CZF_media(m) = mean(CZF);
    CMMSE_media(m) = mean(CMMSE);
    CDPC_media(m) = mean(CDPC);
    
    CMRT_media_iid(m) = mean(CMRT_iid);
    CZF_media_iid(m) = mean(CZF_iid);
    CMMSE_media_iid(m) = mean(CMMSE_iid);
    CDPC_media_iid(m) = mean(CDPC_iid);
    
end
toc

figure
plot(ntx,CMRT_media,'-o','LineWidth',1.5)
hold on
plot(ntx,CZF_media,'r-o','LineWidth',1.5)
plot(ntx,CMMSE_media,'k-o','LineWidth',1.5);
plot(ntx,CDPC_media,'o-','Color',[0,1,0.9],'LineWidth',1.5);
grid on
plot(ntx,CMRT_media_iid,'--b','LineWidth',1.5)
plot(ntx,CZF_media_iid,'--r','LineWidth',1.5)
plot(ntx,CMMSE_media_iid,'--k','LineWidth',1.5)
plot(ntx,CDPC_media_iid,'--','Color',[0,1,0.9],'LineWidth',1.5)
xlabel('Número de antenas en BS (M)')
ylabel('Eficiencia Espectral (bits/s/Hz)');
grid on;
if DL ==1
    legend('MRT','ZF','MMSE','DPC via IWF');
    str = sprintf('Eficiencia espectral frente a M para un enlace descendente K = %d; SNR = %d (dB)',nrx,SNRdB);
else
    legend('MRC','ZF','MMSE','Sum-rate óptimo');
    str = sprintf('Eficiencia espectral frente a M para un enlace ascendente K = %d; SNR = %d (dB)',nrx,SNRdB);
end
title(str);

