function result = GPU_dqdt(h,h2,hx,hZ,hxz,hxx,hzz,q,qx,qz,qxx,qzz,P,pz,px,pxz,delta,eta,zeta,datahLx,corr)
%the derivative of the streamwise flow rate by time by using the simplified model. 
%Use when we use the GPU version of runSimulation
%  

result=  (((9./7).*q.*q./(h.*h) - (5.*cot_theta./2)./re_num.*h).*hx - (18.0./7.0)./h.*q.*qx + q./(7.*h).*qx + ...
                  (5./(2.*delta_num*re_num).*(h - q./h2)) + (5.*delta_num.*delta_num.*weber_num)./(2.*re_num).*h.*real(datahLx) + ...
                  delta_num./re_num.*(9./2.*qxx - 9./2./h.*qx.*hx + 4.*q./h2.*hx.*hx - 6.*q./h.*hxx) + ...
                  (5.*deb_num.*(1 - r_vis))./(2.*re_num).*(6.*q.*qx./h3 - 6.*q.*q.*hx./(h2.*h2)))./(1 - (5.*deb_num.*(1 - r_vis))./(2.*re_num.*h2));
end