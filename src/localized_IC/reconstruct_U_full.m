function [datah, dataq, dataS1, dataS2]   = reconstruct_U_full(U,N1,N2,M1,M2)
%Extracting Data from U Vector
sizeU=size(U(end,:));

% if sizeU(2) == M2*(N1+1)*7
    reconstr=reshape(U(end,:),(N1+1)*M2,7);  %Full second-order model

    datah = reshape(reconstr(:,1),M2,N1+1)';
    dataq = reshape(reconstr(:,2),M2,N1+1)';
    % datap = reshape(reconstr(:,3),M2,N1+1)';
    dataS1 = reshape(reconstr(:,4),M2,N1+1)';
    dataS2 = reshape(reconstr(:,5),M2,N1+1)';
    % dataR1 = reshape(reconstr(:,6),M2,N1+1)';
    % dataR2 = reshape(reconstr(:,7),M2,N1+1)';

    %Reconstruct h,q,p by inverse fourier transformation
    datah   = ffti(datah*M1*M2,N1,N2,M1,M2);
    dataq   = ffti(dataq*M1*M2,N1,N2,M1,M2);
    % datap   = ffti(datap*M1*M2,N1,N2,M1,M2);
    dataS1 = ffti(dataS1*M1*M2,N1,N2,M1,M2);
    dataS2 = ffti(dataS2*M1*M2,N1,N2,M1,M2);
    % dataR1 = ffti(dataR1*M1*M2,N1,N2,M1,M2);
    % dataR2 = ffti(dataR2*M1*M2,N1,N2,M1,M2);
% else
%     reconstr=reshape(U(end,:),(N1+1)*M2,3);  %Simplified model

%     datah=reshape(reconstr(:,1),M2,N1+1)';
%     dataq=reshape(reconstr(:,2),M2,N1+1)';
%     datap=reshape(reconstr(:,3),M2,N1+1)';

%     %Reconstruct h,q,p by inverse fourier transformation
%     datah   = ffti(datah*M1*M2,N1,N2,M1,M2);
%     dataq   = ffti(dataq*M1*M2,N1,N2,M1,M2);
%     datap   = ffti(datap*M1*M2,N1,N2,M1,M2);
% end

% datah=reshape(reconstr(:,1),M2,N1+1)';
% dataq=reshape(reconstr(:,2),M2,N1+1)';
% datap=reshape(reconstr(:,3),M2,N1+1)';

%Reconstruct h,q,p by inverse fourier transformation
% datah   = ffti(datah*M1*M2,N1,N2,M1,M2);
% dataq   = ffti(dataq*M1*M2,N1,N2,M1,M2);
% datap   = ffti(datap*M1*M2,N1,N2,M1,M2);
end

