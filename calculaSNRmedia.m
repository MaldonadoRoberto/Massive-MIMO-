function [SNR_media,SNR_media_dB] = calculaSNRmedia(ch_par,Pt,No)

% Obtenemos el numero de tx y rx;

K = size(ch_par.rx_array,2);
M = ch_par.tx_array.no_elements;


path_loss_dB = ch_par.get_pl;  % Se obtienen los vlaores de path loss en dB

% Pasamos el pathloss a unidades lineales
path_loss = 10.^(path_loss_dB./10);

for i = 1:K
    
    PL = path_loss(i);
    
    P_rec(i) = Pt/PL;
    
    SNR_rec(i) = P_rec(i)/No;
    
    SNR_recDB = 10*log10(SNR_rec(i));
    
     
end

SNR_media = mean(SNR_rec);
SNR_media_dB = mean(SNR_recDB);

