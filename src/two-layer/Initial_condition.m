function [datah,datap,dataq,U] = Initial_condition(stNoise, Lz, Lx, h0, Ax, k0x, Bx, nwx, k0z, Az, nwz, Bz, spNoise, NN1, NN2)
%calculate the initial condition for the full model
%Here comes the parameter change
N1=2^NN1; M1=2*N1; N2=2^NN2; M2=2*N2; N=4*N2*(N1+1);

datah=[];
dataq=[];
datap=[];

%%%%% TL modified here %%%%%
%%%%% problem-specific parameters (dimensionless) %%%%%
%%%%% ρ = 0.5, m = 0.25, h1 = 0.5, Fr = 0.80, So = 0.05 %%%%%
% k0x = parameter.Shkadov_k0x;
So = 0.050;
k0z = parameter.Shkadov_k0z;
% k0x = 2*pi./(8.5./So);
k0x = 2*pi./(100.0);
% h0 = parameter.h0;
h1 = 0.50;
m = 0.25;
rho = 0.50;
fr = 0.80;
we = 0.00;

r = (1.0-h1)/h1;
beta = (m+r).^(-1).*(1+r.*rho);

for n=[1:M2]
    randz(n)=spNoise*rand(1);
end

    for m=[1:M1]
       
       x=(m-1)*Lx/M1;
       for n=[1:M2]
            z=(n-1)*Lz/M2;
            datah(m,n)=h0*(1+Ax*sin(nwx*k0x*x)+Bx*cos(nwx*k0x*x)+Az*sin(nwz*k0z*z)+Bz*cos(nwz*k0z*z)+stNoise*rand(1))+randz(n);
            datap(m,n)=0;
            % power-law consideration here %
            % dataq(m,n)=datah(m,n).^3.*phi0n;
            %%%%% TL modified here %%%%%
            %%%%% we keep the "constant-q" configuration %%%%%
            dataq(m,n)=m.*(m+r.*(4+3.*r.*rho)).*(m.^2+r.^4.*rho+m.*r.*(4+3.*r+r.*(3+4.*r).*rho)).^(-1);
       end
    end

%/*calculate mass*/
mass=sum(sum(datah));
mass=mass/(M1*M2);

%/*set mass to h0 for random function*/
if(spNoise~=0)
    datah=datah-mass+h0;
end

%/*Fourier transform of the initial condition*/
datah=fft2(datah);
datap=fft2(datap);
dataq=fft2(dataq);

%Reciprocity
datah(1,N2+2:M2)=fliplr(conj(datah(1,2:N2)));
datap(1,N2+2:M2)=fliplr(conj(datap(1,2:N2)));
dataq(1,N2+2:M2)=fliplr(conj(dataq(1,2:N2)));

datah=conj(datah);
datap=conj(datap);
dataq=conj(dataq);

datahSV=reshape(datah',1,[]);
datapSV=reshape(datap',1,[]);
dataqSV=reshape(dataq',1,[]);

M1M2=M1*M2;

U=[datahSV(1:N/2) dataqSV(1:N/2) datapSV(1:N/2)]./M1M2;

datah=datah./M1M2;
dataq=dataq./M1M2;
datap=datap./M1M2;
end