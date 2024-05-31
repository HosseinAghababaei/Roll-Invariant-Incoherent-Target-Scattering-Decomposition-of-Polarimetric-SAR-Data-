


clc
clear 
load('T_Flevoland.mat')

[Alphas, Phis,Taw, Ori_max, M] = Aghababaei_decomposion(T);


figure();
subplot(231); imagesc(Alphas);
subplot(232); imagesc(Phis);
subplot(233); imagesc(Taw);
subplot(234); imagesc(Ori_max);
subplot(235); imagesc(M,[0 1]);
