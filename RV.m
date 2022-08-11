%% model comparison 
function [TRV,RV_num,sigma2]=RV(res_1,res_3,Z_ss)
%{ 
% model inputs are residuals from model 1 and model 3, res_1, res_3, and instrument variables Z_ss.  
res_1=pr-1000*downmarkup_n-X_s*gamma_n_1;
res_2=pr-1000*upmarkup_n-1000*downmarkup_n-X_s*gamma_n_2;
res_3=pr-1000*upmarkup_n3-1000*downmarkup_n-X_s*gamma_n_3;
res_4=pr-1000*downmarkup_fullcollusion_n-X_s*gamma_n_4;
res_5=pr-1000*downmarkup_allcollusion-X_s*gamma_n_5;
res_6=pr-1000*dm-X_s*gamma6;
res_666=pr-1000*dm-X_s*gamma_n_666;
%}

%compare 3 and 1
%given res_1 and res_3

n=size(Z_ss,1); %number of the data
L = size(Z_ss,2);%number of the instrument

g=Z_ss'*res_1/n;
h=Z_ss'*res_3/n;

r1=res_1.*Z_ss;
W1=inv(r1'*r1/n);
v_1=r1'*r1/n;gmm1=g'*inv(v_1)*g;

r3=res_3.*Z_ss;
W3=inv(r3'*r3/n);
v_3=r3'*r3/n;gmm3=h'*inv(v_3)*h;

%si 1

si1=sqrtm(W1)*(r1'-g)-1/2*(sqrtm(sqrtm(W1))*sqrtm(sqrtm(W1))*sqrtm(sqrtm(W1)))*...
    ((Z_ss*(sqrtm(sqrtm(W1))*sqrtm(sqrtm(W1))*sqrtm(sqrtm(W1))*g)).*Z_ss)'-...
    1/2*sqrtm(W1)*g;

%si 2

si2=sqrtm(W3)*(r3'-h)-1/2*(sqrtm(sqrtm(W3))*sqrtm(sqrtm(W3))*sqrtm(sqrtm(W3)))*...
    ((Z_ss*(sqrtm(sqrtm(W3))*sqrtm(sqrtm(W3))*sqrtm(sqrtm(W3))*h)).*Z_ss)'-...
    1/2*sqrtm(W3)*h;

V11=1/n*si1*si1';V12=1/n*si1*si2';V22=1/n*si2*si2';

sigma2=4*(g'*sqrtm(W1)*V11*sqrtm(W1)*g+h'*sqrtm(W3)*V22*sqrtm(W3)*h-...
    2*g'*sqrtm(W1)*V12*sqrtm(W3)*h);
RV_num=(n^0.5)*(g'*inv(v_1)*g - h'*inv(v_3)*h);
TRV = RV_num/sqrt(sigma2);




