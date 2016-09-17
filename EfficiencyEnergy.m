%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      Script: EfficiencyEnergy                     %%%%
%%%%                                                                   %%%%
%%%%        Script la relación entre la eficiencia espectral de un     %%%%
%%%%        sistema y la eficiencia energetica definida como la        %%%%
%%%%        eficiencia esptral entre la potencia transmitida.          %%%%
%%%%       Por lo tanto, se obtiene la eficiencia espectral realizando %%%%
%%%%        variaciones en la SNR y para dos configuraciones distantas %%%%
%%%%        de MIMO masivo.                                            %%%%
%%%%                                                                   %%%%
%%%%           En este script se utiliza la normalización de           %%%%
%%%%           de las matrices de canal para que solo contengan        %%%%
%%%%               los efectos de small scale fading                   %%%%
%%%%                                                                   %%%%
%%%%  No se calculan los resultados para DPC (alto tiempo de ejecución)%%%%  
%%%%                                                                   %%%%
%%%%                  Resultados mostrados en 7.4.3                    %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc
addpath('..')

% Parametros
nrx = 10;
precision = 5;
iteraciones = 2;
realizaciones = 25;
SNRdB = [-10:5:20];
M = [50, 100];
DL = 1;

for m = 1:length(M)
    fprintf('M = %d',M(m));
    for n = 1:length(SNRdB)
        
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
            
            % Debido al numero tan alto realizaciones se van a coger 5 al azar
            % que seran las utilizadas para obtener la eficiencia espectral;
            
            index = randperm(realizaciones,10);
            % Calcular el sum rate para cada realización del canal y realizar la
            % media:
            
            for k = 1:length(index)
                
                H = canales_normalizados{index(k)};
                H_iid = canales_normalizados_iid{index(k)};
                
                % Redución en un factor M de la potencia transmitida;

                % SNRdB_reducido = SNRdB;
                if DL == 1
                    [C_mrt(k),C_zf(k),C_mmse(k)] = sum_rate_DL_sinDPC(H,SNRdB(n));
                    [C_mrt_iid(k),C_zf_iid(k),C_mmse_iid(k)] = sum_rate_DL_sinDPC(H_iid,SNRdB(n));
                else
                    [C_mrt(k),C_zf(k),C_mmse(k),C_dpc(k)] = sum_rate_UL(H,SNRdB(n));
                    [C_mrt_iid(k),C_zf_iid(k),C_mmse_iid(k),C_dpc_iid(k)] = sum_rate_UL(H_iid,SNRdB(n));
                end
            end
            
            CMRT(i) = mean(C_mrt);
            CZF(i) = mean(C_zf);
            CMMSE(i) = mean(C_mmse);
            
            CMRT_iid(i) = mean(C_mrt_iid);
            CZF_iid(i) = mean(C_zf_iid);
            CMMSE_iid(i) = mean(C_mmse_iid);
            
            
        end
        
        eval(sprintf('mrt_%d(n) = mean(CMRT);',M(m)));  % Sumatoria del rate de todos los usuarios para cada iteración
        eval(sprintf('zf_%d(n) = mean(CZF);',M(m)));  % Sumatoria del rate de todos los usuarios para cada iteración
        eval(sprintf('mmse_%d(n) = mean(CMMSE);',M(m)));  % Sumatoria del rate de todos los usuarios para cada iteración
        eval(sprintf('mrt_iid_%d(n) = mean(CMRT_iid);',M(m)));  % Sumatoria del rate de todos los usuarios para cada iteración
        eval(sprintf('zf_iid_%d(n) = mean(CZF_iid);',M(m)));  % Sumatoria del rate de todos los usuarios para cada iteración
        eval(sprintf('mmse_iid_%d(n) = mean(CMMSE_iid);',M(m)));  % Sumatoria del rate de todos los usuarios para cada iteración
        fprintf('nueva SNR %d\n',SNRdB(n))
    end
    
    
end

% Eficiencia energetica:
EE_50 = mrt_50./(10.^(SNRdB./10))
EE_100 = mrt_100./(10.^(SNRdB./10))
EE_zf_50 = zf_50./(10.^(SNRdB./10))
EE_zf_100 = zf_100./(10.^(SNRdB./10))
EE_mmse_50 = mmse_50./(10.^(SNRdB./10))
EE_mmse_100 = mmse_100./(10.^(SNRdB./10))

EE_iid_50 = mrt_iid_50./10.^(SNRdB./10)
EE_iid_100 = mrt_iid_100./10.^(SNRdB./10)
EE_iid_zf_50 = zf_iid_50./10.^(SNRdB./10)
EE_iid_zf_100 = zf_iid_100./10.^(SNRdB./10)
EE_iid_mmse_50 = mmse_iid_50./10.^(SNRdB./10)
EE_iid_mmse_100 = mmse_iid_100./10.^(SNRdB./10)


figure
semilogy(mrt_50,EE_50,'o--','LineWidth',1)
hold on
semilogy(zf_50,EE_zf_50,'o--k','LineWidth',1)
semilogy(mmse_50,EE_mmse_50,'o--r','LineWidth',1)
semilogy(mrt_100,EE_100,'o-','LineWidth',1.5)
semilogy(zf_100,EE_zf_100,'o--k','LineWidth',1.5)
semilogy(mmse_100,EE_mmse_100,'o--r','LineWidth',1.5)
grid on
xlabel('Eficiencia espectral (bits/s/Hz)')
ylabel('Eficiencia energética (bits/J/Hz)')
if DL ==1
    legend('MRT','ZF','MMSE');
    str = sprintf('Eficiencia espectral vs eficiencia energética para M = 50 y M = 100 (DL)',nrx,SNRdB);
else
    legend('MRC','ZF','MMSE');
    str = sprintf('Eficiencia espectral vs eficiencia energética para M = 50 y M = 100 (UL)',nrx,SNRdB);
end
title(str);

