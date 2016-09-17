function [norm_channels] = realizaNormalizacion(H)

realizaciones = length(H);

for i = 1:realizaciones
    realizacion = H{i};
    % obtrención de la matriz que contiene la potencia de cad ausuario;
    pow(:,:,i) = abs(realizacion).^2;
end

pow_media = mean(pow,3);
pow_media_media = mean(pow_media);
D = diag(pow_media_media);

for k = 1:realizaciones
    norm_channels{k} = H{k}*D^(-1/2);
end
end