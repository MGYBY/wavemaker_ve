function [storeData] = save_csv(parameter,U,M1,M2,N1,N2,NN2,h0,storeData,kx,kz,ttime,length)
%make all the output for the user

[datah,dataq,~] = reconstruct_U(U,N1,N2,M1,M2);
% hmax=real(max(max(datah)));
% hmin=real(min(min(datah)));
% q = mean(mean(real(dataq)*parameter.Reynolds*3));
% storeData.hmax=[storeData.hmax hmax];
% storeData.hmin=[storeData.hmin hmin];
% storeData.time=[storeData.time ttime];
% storeData.q=[storeData.q q];
storeData.hmax=[storeData.hmax 0.0];
storeData.hmin=[storeData.hmin 0.0];
storeData.time=[storeData.time 0.0];
storeData.q=[storeData.q 0.0];

% combined_mat = [0:1/(M1-1):1)*ScaleX, flipud(real(datah(:,1)))*ScaleH, flipud(real(dataq(:,1)))*ScaleQ];
ScaleX = 1.0*length;
ScaleH = 1.0;
ScaleQ = 1.0;
ScaleT = 1.0;
mod_x_coord = transpose(0:1.0*ScaleX/(M1-1):1.0*ScaleX);
combined_mat = [mod_x_coord, flipud(real(datah(:,1)))*ScaleH, flipud(real(dataq(:,1)))*ScaleQ];
filename_text = ([parameter.Name '/' sprintf('%g', ttime*ScaleT) 'text_field.csv']);
csvwrite (filename_text, combined_mat);
end

