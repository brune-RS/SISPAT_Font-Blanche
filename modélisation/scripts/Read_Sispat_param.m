%% Read_Sispat_param.m
% Cette fonction lit un fichier de paramètres SiSPAT et le stocke dans une
% structure

function Param_SiSPAT=Read_Sispat_param(file_param)

fid=fopen(file_param,'r+');

lines=textscan(fid,'%s','delimiter','\n');
Param_SiSPAT=lines{1};
fclose(fid)
end