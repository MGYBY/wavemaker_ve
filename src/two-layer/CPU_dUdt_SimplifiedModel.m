function [dUdt] = CPU_dUdt_SimplifiedModel(yU,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N,A,hs)
%function calculating the equation for the simplified model that will be solved by ode45
%
%see also fft, ode45
% 

%reconstruct data from U-vector (filmthickness h, streamwise flow rate q
%and spanwise flow rate p

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

% other important parameters %
re = (fr.^2.0./So)./(((r.*(1+r).*(1+r.*rho))./(m+r)+(m+r.^3.0.*rho)./(3.0.*m))./(4.*(1+r).^3.0)); % the "R" parameter %
G_param = 4.*(1+r).^3.*(r.*(1+r).*(m+r).^(-1).*(1+r.*rho)+(1/3).*m.^(-1).*(m+r.^3.*rho)).^(-1); % the "G" parameter %

epsilon = So;

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

      %/* x-derivatives */
datahx=kx.*imag(datah)-kx.*real(datah)*1i;      
dataqx=kx.*imag(dataq)-kx.*real(dataq)*1i;
datapx=kx.*imag(datap)-kx.*real(datap)*1i;

%Transform variables back to real space
hx  = real(ffti(datahx*M1M2,N1,N2,M1,M2));
px  = real(ffti(datapx*M1M2,N1,N2,M1,M2));
qx  = real(ffti(dataqx*M1M2,N1,N2,M1,M2));

%%%%% TL modified here %%%%%
%%%%% could be simpler in 2D (dummy) %%%%%
      %/* z-derivatives */
% datahz=kz.*imag(datah)-kz.*real(datah)*1i;
% dataqz=kz.*imag(dataq)-kz.*real(dataq)*1i;
% datapz=kz.*imag(datap)-kz.*real(datap)*1i;

% %Transform variables back to real space
% hZ  = real(ffti(datahz*M1M2,N1,N2,M1,M2));
% pz  = real(ffti(datapz*M1M2,N1,N2,M1,M2));
% qz  = real(ffti(dataqz*M1M2,N1,N2,M1,M2));

      %/* Double x-derivatives */
datahxx=-kx2.*real(datah)-kx2.*imag(datah)*1i;
dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*1i;
datapxx=-kx2.*real(datap)-kx2.*imag(datap)*1i;  

%Transform variables back to real space      
hxx = real(ffti(datahxx*M1M2,N1,N2,M1,M2));
pxx = real(ffti(datapxx*M1M2,N1,N2,M1,M2));
qxx = real(ffti(dataqxx*M1M2,N1,N2,M1,M2));
   
%%%%% TL modified here %%%%%
%%%%% could be simpler in 2D (dummy) %%%%%
      %/* Double z-derivatives */
% datahzz=-kz2.*real(datah)-kz2.*imag(datah)*1i;
% dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*1i;
% datapzz=-kz2.*real(datap)-kz2.*imag(datap)*1i; 

% %Transform variables back to real space
% hzz = real(ffti(datahzz*M1M2,N1,N2,M1,M2));
% pzz = real(ffti(datapzz*M1M2,N1,N2,M1,M2));
% qzz = real(ffti(dataqzz*M1M2,N1,N2,M1,M2));


%%%%% TL modified here %%%%%
%%%%% could be simpler in 2D (dummy) %%%%%
     %/* Cross derivatives */
% datahxz=-kxz.*real(datah)-kxz.*imag(datah)*1i;
% dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*1i;
% datapxz=-kxz.*real(datap)-kxz.*imag(datap)*1i;

% %Transform variables back to real space
% hxz = real(ffti(datahxz*M1M2,N1,N2,M1,M2));
% pxz = real(ffti(datapxz*M1M2,N1,N2,M1,M2));
% qxz = real(ffti(dataqxz*M1M2,N1,N2,M1,M2));

    %/* Gradient of the laplacien of h (=curvature K) 
datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*1i);
datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*1i); 
    
%%%%% TL modified here %%%%%
%%%%% could be simpler in 2D %%%%%
% dhdt=-real(dataqx)-real(datapz)-imag(dataqx)*1i-imag(datapz)*1i;
dhdt=-real(dataqx)-imag(dataqx)*1i;

%Transform variables back to real space
h   = real(ffti(datah*M1M2,N1,N2,M1,M2));
P   = real(ffti(datap*M1M2,N1,N2,M1,M2));
q   = real(ffti(dataq*M1M2,N1,N2,M1,M2));


datahLx = ffti(datahLx*M1M2,N1,N2,M1,M2);
datahLz = ffti(datahLz*M1M2,N1,N2,M1,M2);

h2=h.*h;
h3=h2.*h;

%%%%% TL modified here %%%%%
%%%%% dummy %%%%%
% Ineq=0; %Switch to 1 for second order corrections

%%%%% TL modified here %%%%%
a1 = ((-1)+h).^(-1).*h.^(-1).*(1+((-1)+m).*h).^(-1).*(9.*m.*h.^2+3.*((-1)+h.*(2+(-4).*m+((-1)+m).*h)).*q);

b1 = ((-1)+h).^(-1).*h.^(-1).*(1+((-1)+m).*h).^(-1).*(9.*q+(-6).*h.*(2+q)+3.*h.^2.*(4+(-1).*m+((-1)+m).*q));

a0 = 6.*m.*((-1)+h).^(-1).*(1+((-1)+m).*h).^(-1).*((-1).*h+q);

% b0 = (h+((-1)+m).*h.^2).^(-1).*((-6).*h+6.*q);

a1q = 3.*h.^(-1)+(-9).*m.*((-1)+h).^(-1).*(1+((-1)+m).*h).^(-1);

a0q = 6.*m.*((-1)+h).^(-1).*(1+((-1)+m).*h).^(-1);

b1q = 3.*((-1)+h).^(-1)+(-9).*h.^(-1)+9.*((-1)+m).*(1+((-1)+m).*h).^(-1);

b0q = 6.*(h+((-1)+m).*h.^2).^(-1);

a0x = 6.*m.*((-1)+h).^(-2).*(1+((-1)+m).*h).^(-2).*((1+((-1)+m).*h.^2+((-2)+m).*q+(-2).*((-1)+m).*h.*q).*hx+((-1)+h).*(1+((-1)+m).*h).*qx);

a1x = ((-1)+h).^(-2).*h.^(-2).*(1+((-1)+m).*h).^(-2).*((-3).*(3.*m.*h.^2.*(1+((-1)+m).*h.^2)+(1+h.*(2.*((-2)+m)+h.*(6+4.*((-3)+m).*m+((-1)+m).*h.*(4+(-8).*m+((-1)+m).*h)))).*q).*hx+3.*((-1)+h).*h.*(1+((-1)+m).*h).*((-1)+h.*(2+(-4).*m+((-1)+m).*h)).*qx);

b0x = (h+((-1)+m).*h.^2).^(-2).*(6.*(((-1)+m).*h.*(h+(-2).*q)+(-1).*q).*hx+6.*h.*(1+((-1)+m).*h).*qx);

b1x = 3.*((-1)+h).^(-2).*h.^(-2).*(1+((-1)+m).*h).^(-2).*(3.*q.*hx+3.*h.*(2.*((-2)+m).*q.*hx+(-1).*qx)+((-1)+m).^2.*h.^5.*qx+h.^2.*((4+(-3).*m+2.*(7+(-6).*m).*q).*hx+(8+(-3).*m).*qx)+((-1)+m).*h.^4.*(((-4)+m+q+(-1).*m.*q).*hx+(-1).*m.*qx)+2.*h.^3.*(2.*((-1)+m).*(2+q).*hx+((-3)+2.*m).*qx));

a0xx = 6.*m.*((-1)+h).^(-3).*(1+((-1)+m).*h).^(-3).*((-2).*(2+(-1).*m+((-1)+m).^2.*h.^3+(-1).*(3+((-3)+m).*m).*q+(-3).*((-1)+m).^2.*h.^2.*q+3.*((-1)+m).*h.*(1+((-2)+m).*q)).*hx.^2+(-2).*((-1)+h).*(1+((-1)+m).*h).*(2+(-1).*m+2.*((-1)+m).*h).*hx.*qx+((-1)+h).*(1+((-1)+m).*h).*((1+((-1)+m).*h.^2+((-2)+m).*q+(-2).*((-1)+m).*h.*q).*hxx+((-1)+h).*(1+((-1)+m).*h).*qxx));

a1xx = 3.*((-1)+h).^(-3).*h.^(-3).*(1+((-1)+m).*h).^(-3).*((-2).*q.*hx.^2+h.*(2.*hx.*qx+q.*((-6).*((-2)+m).*hx.^2+hxx))+((-1)+m).^3.*h.^8.*qxx+(-1).*h.^2.*((-6).*((-2)+m).*hx.*qx+3.*q.*(2.*(5+((-5)+m).*m).*hx.^2+(-1).*((-2)+m).*hxx)+qxx)+h.^3.*(((-6).*((-2)+m).*m+2.*(20+m.*((-39)+(21+(-4).*m).*m)).*q).*hx.^2+6.*((-1)+m).*((-5)+2.*m).*hx.*qx+3.*(m+((-1)+m).*((-5)+2.*m).*q).*hxx+(-6).*((-1)+m).*qxx)+(-1).*((-1)+m).^2.*h.^7.*(2.*((-1)+m).*hx.*qx+(3.*m+((-1)+m).*q).*hxx+6.*((-1)+m).*qxx)+h.^4.*(6.*((-1)+m).*(3.*m+(5+m.*((-11)+4.*m)).*q).*hx.^2+4.*((-10)+m.*(24+m.*((-15)+2.*m))).*hx.*qx+(3.*((-2)+m).*m+2.*((-10)+m.*(24+m.*((-15)+2.*m))).*q).*hxx+(-3).*(5+3.*((-3)+m).*m).*qxx)+((-1)+m).*h.^6.*(2.*((-1)+m).*(3.*m+((-1)+m).*q).*hx.^2+6.*((-1)+m).*((-2)+3.*m).*hx.*qx+3.*(((-2)+m).*m+((-1)+m).*((-2)+3.*m).*q).*hxx+3.*(5+3.*((-3)+m).*m).*qxx)+h.^5.*(3.*((-1)+m).*((-2).*(5+m.*((-11)+4.*m)).*hx.*qx+q.*((-4).*((-1)+m).*((-1)+2.*m).*hx.^2+((-5)+(11+(-4).*m).*m).*hxx))+2.*(10+m.*((-24)+(15+(-2).*m).*m)).*qxx));

b0xx = (h+((-1)+m).*h.^2).^(-3).*(2.*(1+2.*((-1)+m).*h).^2.*((-6).*h+6.*q).*hx.^2+(-2).*(1+2.*((-1)+m).*h).*(h+((-1)+m).*h.^2).*hx.*((-6).*hx+6.*qx)+(-1).*(h+((-1)+m).*h.^2).*((-6).*h+6.*q).*(2.*((-1)+m).*hx.^2+(1+2.*((-1)+m).*h).*hxx)+(h+((-1)+m).*h.^2).^2.*((-6).*hxx+6.*qxx));

b1xx = 3.*((-1)+h).^(-3).*h.^(-3).*(1+((-1)+m).*h).^(-3).*(6.*q.*hx.^2+3.*h.*((-2).*hx.*qx+q.*(6.*((-2)+m).*hx.^2+(-1).*hxx))+((-1)+m).^3.*h.^8.*qxx+3.*h.^2.*((-6).*((-2)+m).*hx.*qx+3.*q.*(2.*(5+((-5)+m).*m).*hx.^2+(-1).*((-2)+m).*hxx)+qxx)+(-1).*h.^3.*(2.*(4+3.*((-2)+m).*m+(56+27.*((-3)+m).*m).*q).*hx.^2+2.*(41+(-39).*m+6.*m.^2).*hx.*qx+(4+(-3).*m+(41+(-39).*m+6.*m.^2).*q).*hxx+2.*(7+(-3).*m).*qxx)+(-1).*((-1)+m).^2.*h.^7.*(2.*((-1)+m).*hx.*qx+(4+(-1).*m+((-1)+m).*q).*hxx+2.*((-1)+m).*qxx)+h.^4.*(6.*((-1)+m).*((-4)+3.*m+((-11)+9.*m).*q).*hx.^2+4.*(22+(-30).*m+9.*m.^2).*hx.*qx+(16+3.*((-6)+m).*m+2.*(22+(-30).*m+9.*m.^2).*q).*hxx+(25+3.*((-7)+m).*m).*qxx)+h.^5.*(3.*((-1)+m).*((-4).*((-1)+m).*(2+q).*hx.^2+2.*(7+(-5).*m).*hx.*qx+(8+(-4).*m+(7+(-5).*m).*q).*hxx)+(-2).*(10+3.*((-4)+m).*m).*qxx)+((-1)+m).*h.^6.*(2.*((-1)+m).*(4+(-1).*m+((-1)+m).*q).*hx.^2+2.*((-2)+m+m.^2).*hx.*qx+((-16)+(-1).*((-14)+m).*m+((-2)+m+m.^2).*q).*hxx+((-5)+m+m.^2).*qxx));

Q0 = ((h-1).*(h.^2.*G_param.*(rho-1) + a1)-(m.*b1.*h.^2)./(h-1)).*(h.*(m-1)+1);

Q1 = (1./20)*(h.^2.0).*(8.*h.^2.*(m-1) + 16.*h - 7.*m.*h - 8).*a1q + (1./8).*(h.^2.0).*(5.*h.^2.*(m-1) - 4.*m.*h + 10.*h - 5).*a0q - (1./8).*h.*rho.*(h-1).*(5.*h.^2.*(m-1) + 6.*h - 1).*b0q - (1./20).*h.*rho.*(h-1).*(8.*h.^2.*(m-1) + 9.*h -1).*b1q;

Q1_tilde = (1./re).*h.^2.*(1-h).*(h.*(m-1)+1).*((rho-1).*G_param.*(1./So).*hx + epsilon.^2.*we.*real(datahLx)) + (1./120).*h.*a0.*(h.*a1x-2.*hx.*a1).*(13.*h.^2.*(m-1)-10.*m.*h+26.*h-13) - (3./(40.*h.*m.^2)).*rho.*(h-1).*(3.*h.^2.*(m-1)+4.*h-1).*(m.*h.*b0x.*(m.*h.*b1+a0.*(h-1))-hx.*a0.^2.*(h-1)) - (1./840).*h.*hx.*(a1.^2).*(33.*h.^2.*(m-1)-34.*m.*h+66.*h-33) + (1./840).*rho.*h.*hx.*(b1.^2).*(33.*h.^2.*(m-1)+32.*h+1) + (3./40).*h.*(h.*a0x.*a0+h.*a0x.*a1-hx.*a0.^2).*(3.*h.^2.*(m-1)-2.*m.*h+6.*h-3) + (1./(60.*m)).*rho.*(h-1).*a0.*(13.*h.^2.*(m-1)+16.*h-3).*(hx.*b1-(1./2).*(h-1).*b1x) - (1./840).*rho.*h.*(h-1).*(111.*h.^2.*b1.*(m-1)+142.*b1h-31.*b1).*b1x + (1./840).*h.^2.*(a1.*a1x).*(111.*h.^2.*(m-1)+222.*h-80.*m.*h-111) + qx.*((-1./40).*h.*rho.*(7.*h.^2.*(m-1)+6.*h+1).*b1+(1./20).*h.*rho.*(h-1).*b1h.*(8.*h.^2.*(m-1)+9.*h-1)-(1./20).*h.^2.*a1h.*(8.*h.^2.*(m-1)+16.*h-7.*m.*h-8)+(1./40).*h.*a1.*(7.*h.^2.*(m-1)-8.*m.*h+14.*h-7)+(1./(8.*m).*rho.*(h-1).*(5.*h.^2.*(m-1)+6.*h-1).*(m.*h.*b0h-a0))-(1./8).*h.*(5.*h.^2.*(m-1)-4.*m.*h+10.*h-5).*(a0h.*h-a0));

Q2 = (3.*h.^2.*(m-1)+4.*h-1).*(3./5.*h.*m.*b1xx.*(h-1)+m.*h.*b0xx.*(h-1)-1./5.*m.*h.*b1.*hxx-a0.*h.*hxx+a0.*hxx) - (1./5.*h.*(3.*h.^2.*(m-1)-2.*m.*h+6.*h-3).*(3.*h.*a1xx+5.*h.*a0xx-5.*a0.*hxx-a1.*hxx)) - hx.*(h.^2.*(m-1)+3.*m.*h+2.*h-1).*(a0x.*h-a0.*hx) + (1./5).*(h.^2.*(m-1)+2.*h+m.*h-1).*a1.*(hx.^2) + hx.*(h.^2.*(m-1)+5.*h-4).*(b0x.*m.*h-a0.*hx) - (1./5.*h.*hx.*a1x).*(9.*h.*(m-1)+4.*m.*h+18.*h-9) + (1./5.*h.*m./(1-h).*b1.*(hx.^2).*(h.^2.*(m-1)-2-h.^2)) + (1./5.*h.*m.*hx.*b1x.*(9.*h.^2.*(m-1)+22.*h-13));

dqdt = -1.0.*Q0./(epsilon.*re.*Q1) - Q1_tilde./Q1 - epsilon.*Q2./(re.*Q1);
% or first-order model
% dqdt = -1.0.*Q0./(epsilon.*re.*Q1) - Q1_tilde./Q1;

%%%%% TL modified here %%%%%
%%%%% could be simpler in 2D %%%%%
dpdt = h.*0.00;

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