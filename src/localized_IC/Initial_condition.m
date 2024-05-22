function [datah,datap,dataq,U] = Initial_condition(stNoise, Lz, Lx, h0, Ax, k0x, Bx, nwx, k0z, Az, nwz, Bz, spNoise, NN1, NN2)
%calculate the initial condition for the full model
%Here comes the parameter change
N1=2^NN1; M1=2*N1; N2=2^NN2; M2=2*N2; N=4*N2*(N1+1);

datah=[];
dataq=[];
datap=[];

for n=[1:M2]
    randz(n)=spNoise*rand(1);
end

    for m=[1:M1]
       
       x=(m-1)*Lx/M1;
       for n=[1:M2]
            z=(n-1)*Lz/M2;
            % account for localized dist
            if (x<=Lx/nwx && x>=Lx/nwx/2.0)
                % note that k0x is the wavenumber for the localized dist
                datah(m,n)=h0*(1+Ax*sin(2.0*pi*(x-Lx/nwx/2.0)/(Lx/nwx))+Bx*cos(nwx*k0x*x)+Az*sin(nwz*k0z*z)+Bz*cos(nwz*k0z*z)+stNoise*rand(1))+randz(n);
                datap(m,n)=0;
                dataq(m,n)=datah(m,n).^3/3;
                % dataq(m,n)=(((3.*datah(m,n).^2)/2-datah(m,n)^3./2)./datah(m,n)).*datah(m,n)
            else
                datah(m,n)=h0*1.0;
                datap(m,n)=0;
                dataq(m,n)=datah(m,n).^3/3;
            end
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
