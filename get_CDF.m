%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      Script: getCDF                               %%%%
%%%%                                                                   %%%%
%%%%        Script en el que se obtienen las CDFs de un vector x       %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CDF,intervalos] = get_CDF(x,int)

[PDF,intervalos] = hist(x,int); % Generación del histograma
% dividido en int intervalos;

PDF = PDF/size(x,2);  % Se normaliza la PDF

for i = 1:int
    CDF(1,i) = sum(PDF([1:i])); % Probabilidad acumulada
end

end

