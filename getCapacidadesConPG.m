%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      Script: getCapacidadesConPG                  %%%%
%%%%                                                                   %%%%
%%%%        Script en el que se obtiene la eficiencia espectral de     %%%%
%%%%        una configuración MIMO en la que en este caso, se tiene    %%%%
%%%%        en cuenta el path gain de los usuarios.                    %%%%
%%%%                                                                   %%%%
%%%%           En este script se utiliza la normalización de           %%%%
%%%%           de las matrices de canal para que solo contengan        %%%%
%%%%               los efectos de small scale fading                   %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%                  Resultados mostrados en 7.4.1                    %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

precision = 5;
ntx=100;
nrx = 8;
Niter = 25;
realizaciones = 20;
Pt = 40;
No = 10^(-100/10)*1e-3; %(100dbm);

for i=1:Niter
   
    fprintf('Iteracion %d de %d\n',i,Niter);
    [H_struct,H_iid_struct,ch_par] = generaCanalQuadriga(precision,ntx,nrx,realizaciones,1,0,0);
    [SNR_media(i), SNR_mediaDB(i)] = calculaSNRmedia(ch_par,Pt,No);
    SNR_text = sprintf('SNR media es de %.2f (%.2f dB)',SNR_media(i),SNR_mediaDB(i));
    disp(SNR_text);
    
    index = randperm(realizaciones,10);
    for j=1:length(index)
        [C_mrt(j,i),C_zf(j,i),C_mmse(j,i),C_dpc(j,i)] = sum_rate_DL_sin_norm(H_struct{index(j)},Pt,No);
        [C_mrt_iid(j,i),C_zf_iid(j,i),C_mmse_iid(j,i),C_dpc_iid(j,i)] = sum_rate_DL_sin_norm(H_iid_struct{index(j)},Pt,No);
    end
    
end

CMRT = mean(C_mrt)
CZF = mean(C_zf);
CMMSE = mean(C_mmse);
CDPC = mean(C_dpc);

CMRT_iid = mean(C_mrt_iid)
CZF_iid = mean(C_zf_iid);
CMMSE_iid = mean(C_mmse_iid);
CDPC_iid = mean(C_dpc_iid);


%% Obtenemos las CDF para poder representar los valores obtenidosy realizar un analisis de los mismos:
[CDF_mrt,int_mrt] = get_CDF(CMRT,50);
[CDF_zf,int_zf] = get_CDF(CZF,50);
[CDF_mmse,int_mmse] = get_CDF(CMMSE,50);
[CDF_DPC,int_DPC] = get_CDF(CDPC,50);

[CDF_mrt_iid,int_mrt_iid,] = get_CDF(CMRT_iid,50);
[CDF_zf_iid,int_zf_iid] = get_CDF(CZF_iid,50);
[CDF_mmse_iid,int_mmse_iid] = get_CDF(CMMSE_iid,50);
[CDF_DPC_iid,int_DPC_iid] = get_CDF(CDPC_iid,50);

p95MRC = prctile(CMRT,95);
p95ZF = prctile(CZF,95);
p95MMSE = prctile(CMMSE,95);
p95DPC = prctile(CDPC,95);

figure
plot(int_mrt,CDF_mrt,'LineWidth',1.5)
hold on
plot(int_zf,CDF_zf,'r','LineWidth',1.5)
plot(int_mmse,CDF_mmse,'k','LineWidth',1.5)
plot(int_DPC,CDF_DPC,'Color',[0,1,0.9],'LineWidth', 1.5);
plot(int_mrt_iid,CDF_mrt_iid,'--b','LineWidth',1.5)
plot(int_zf_iid,CDF_zf_iid,'--r','LineWidth',1.5)
plot(int_mmse_iid,CDF_mmse_iid,'--k','LineWidth',1.5)
plot(int_DPC_iid,CDF_DPC_iid,'--','Color',[0,1,0.9],'LineWidth',1.5);
plot(p95MRC,0.95,'s','Color',[1,0.81,0.1])
plot(p95ZF,0.95,'s','Color',[1,0.81,0.1])
plot(p95MMSE,0.95,'s','Color',[1,0.81,0.1])
plot(p95DPC,0.95,'s','Color',[1,0.81,0.1])
axis([min(min(int_mrt))-5 max(max(int_DPC_iid))+5 0 1.02]);
grid on;
legend('MRT','ZF','MMSE','DPC')
xlabel('Eficiencia Espectral (bits/s/Hz)');
ylabel('P(SP<x)')
str = sprintf('CDF de la capacidad en DL con %d TX antenas y %d usuarios. SNR media = %.2f dB',ntx,nrx,mean(SNR_mediaDB));
title(str);



