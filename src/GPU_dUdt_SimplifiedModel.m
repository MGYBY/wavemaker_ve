function [dUdt] = CPU_dUdt_SimplifiedModel(yU,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N,A,hs)
%function calculating the equation for the simplified model that will be solved by ode45
%
%see also fft, ode45
% 

%reconstruct data from U-vector (filmthickness h, streamwise flow rate q
%and spanwise flow rate p

r_vis = 1.00;
deb_num = 0.0;
re_num = 5.0;
% cot_theta = 1.0/0.11217;
cot_theta = 1.0/1.0;
froude_num = (1.0/cot_theta*re_num)^0.5;
% weber_num = 20.0*cot_theta*(froude_num*froude_num);
weber_num = 20.0;
delta_num = 1.0/cot_theta;

reconstru=reshape(yU,(N1+1)*M2,3);
datah = reconstru(:,1);
dataq = reconstru(:,2);
datap = reconstru(:,3);
clear reconstru

datah=reshape(datah,M2,N1+1)';
dataq=reshape(dataq,M2,N1+1)';
datap=reshape(datap,M2,N1+1)';

%fix the value of the average film thickness
datah(1,1)=h0+imag(datah(1,1)); %Set average film thickness to fixed value

%Wavenumber values needed for derivatives in Fourier space
kx2=kx.*kx;
kz2=kz.*kz;
kxz=kx.*kz;
k2=kx2+kz2;

M1M2=M1*M2;

datah=gpuArray([datah; [(flipud(datah(2:N1,1))) (rot90(datah(2:N1,2:M2),2))]]);
dataq=gpuArray([dataq; [(flipud(dataq(2:N1,1))) (rot90(dataq(2:N1,2:M2),2))]]);
datap=gpuArray([datap; [(flipud(datap(2:N1,1))) (rot90(datap(2:N1,2:M2),2))]]);

      %/* x-derivatives */
datahx=kx.*imag(datah)-kx.*real(datah)*1i;      
dataqx=kx.*imag(dataq)-kx.*real(dataq)*1i;
datapx=kx.*imag(datap)-kx.*real(datap)*1i;

%Transform variables back to real space
hx  = real(fftiGPU(datahx*M1M2,N1,N2,M1,M2));
px  = real(fftiGPU(datapx*M1M2,N1,N2,M1,M2));
qx  = real(fftiGPU(dataqx*M1M2,N1,N2,M1,M2));

      
      %/* z-derivatives */
datahz=kz.*imag(datah)-kz.*real(datah)*1i;
dataqz=kz.*imag(dataq)-kz.*real(dataq)*1i;
datapz=kz.*imag(datap)-kz.*real(datap)*1i;

%Transform variables back to real space
hZ  = real(fftiGPU(datahz*M1M2,N1,N2,M1,M2));
pz  = real(fftiGPU(datapz*M1M2,N1,N2,M1,M2));
qz  = real(fftiGPU(dataqz*M1M2,N1,N2,M1,M2));

      
      %/* Double x-derivatives */
datahxx=-kx2.*real(datah)-kx2.*imag(datah)*1i;
dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*1i;
datapxx=-kx2.*real(datap)-kx2.*imag(datap)*1i;  

%Transform variables back to real space      
hxx = real(fftiGPU(datahxx*M1M2,N1,N2,M1,M2));
pxx = real(fftiGPU(datapxx*M1M2,N1,N2,M1,M2));
qxx = real(fftiGPU(dataqxx*M1M2,N1,N2,M1,M2));
   
      %/* Double z-derivatives */
datahzz=-kz2.*real(datah)-kz2.*imag(datah)*1i;
dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*1i;
datapzz=-kz2.*real(datap)-kz2.*imag(datap)*1i; 

%Transform variables back to real space
hzz = real(fftiGPU(datahzz*M1M2,N1,N2,M1,M2));
pzz = real(fftiGPU(datapzz*M1M2,N1,N2,M1,M2));
qzz = real(fftiGPU(dataqzz*M1M2,N1,N2,M1,M2));

      
     %/* Cross derivatives */
datahxz=-kxz.*real(datah)-kxz.*imag(datah)*1i;
dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*1i;
datapxz=-kxz.*real(datap)-kxz.*imag(datap)*1i;

%Transform variables back to real space
hxz = real(ffti(datahxz*M1M2,N1,N2,M1,M2));
pxz = real(ffti(datapxz*M1M2,N1,N2,M1,M2));
qxz = real(ffti(dataqxz*M1M2,N1,N2,M1,M2));

    %/* Gradient of the laplacien of h (=curvature K) 
datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*1i);
datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*1i); 
    
dhdt=arrayfun(@GPU_dhdt,dataqx,datapz);

%Transform variables back to real space
h   = real(fftiGPU(datah*M1M2,N1,N2,M1,M2));
P   = real(fftiGPU(datap*M1M2,N1,N2,M1,M2));
q   = real(fftiGPU(dataq*M1M2,N1,N2,M1,M2));


datahLx = fftiGPU(datahLx*M1M2,N1,N2,M1,M2);
datahLz = fftiGPU(datahLz*M1M2,N1,N2,M1,M2);

h2=h.*h;
h3=h2.*h;

Ineq=0; %Switch to 1 for second order corrections
       
dqdt=arrayfun(@GPU_dqdt,h,h2,hx,hZ,hxz,hxx,hzz,q,qx,qz,qxx,qzz,P,pz,px,pxz,delta,eta,zeta,datahLx,corr);
dpdt=arrayfun(@GPU_dpdt,h,h2,h3,hx,hZ,hxz,hxx,hzz,q,qx,qz,qxz,P,pz,px,pxx,pzz,datahLz,delta,eta,zeta);
                          
% /* FOURIER SPACE (attention /(M1*M2)!  */
dqdt=fft2(dqdt);
dpdt=fft2(dpdt);

%    /* Alaising filter removes high frequency modes */
dqdt=dqdt.*k2cut/M1M2;
dpdt=dpdt.*k2cut/M1M2;

%Transform data back to state vector
dhdt=conj(dhdt);
dqdt=conj(dqdt);
dpdt=conj(dpdt);

dhdt=reshape(dhdt',1,[]);
dqdt=reshape(dqdt',1,[]);
dpdt=reshape(dpdt',1,[]);

%Required derivative 
dUdt=[dhdt(1:N/2) dqdt(1:N/2) dpdt(1:N/2)]';

end