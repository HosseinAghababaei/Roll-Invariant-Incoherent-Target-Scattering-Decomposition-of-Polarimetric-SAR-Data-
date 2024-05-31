function [Alphas, Phis,Taw, Ori_max, M] = Aghababaei_decomposion(T)

[~,N,Lazi,Lrng] = size(T);

Alphas = zeros(Lazi*Lrng,1);
Phis = zeros(Lazi*Lrng,1);
Taw = zeros(Lazi*Lrng,1);
Ori_max = zeros(Lazi*Lrng,1);
M = zeros(Lazi*Lrng,1);

[ii,jj] = meshgrid(1:Lazi,1:Lrng);
ii = ii(:);   
jj = jj(:);
dispN = round(Lazi*Lrng/10); 

parfor kk=1:Lazi*Lrng

    [V ,D]=eig(T(:,:,ii(kk), jj(kk)));
    [D ,indx] = sort(diag(D),'descend');
    V=V(:,indx);
    

    % For dominant scatterer
    S =zeros(4,1);
    S(1) = sqrt(D(1))*(V(1,1)+V(2,1))./sqrt(2);
    S(2) = sqrt(D(1))*V(3,1)/sqrt(2);
    S(3) = sqrt(D(1))*V(3,1)/sqrt(2);
    S(4) = sqrt(D(1))*(V(1,1)-V(2,1))/sqrt(2);

    
    alpha=(1/sqrt(2))*(S(1)+S(4));
    ro=S(1)-(alpha/sqrt(2))+ 1i*S(2); % 0.5*(S(1)-S(4)+2*i*S(2))
    lambda= ro-2*1i*S(2); %0.5*(S(1)-S(4)-2*i*S(2))
    a =norm(S);
    cos_taw=sqrt( ( abs(alpha)^2+( (abs(lambda)+abs(ro))^2 )/2 )/(abs(alpha)^2+abs(lambda)^2+abs(ro)^2));
    sin_taw=((abs(ro)-abs(lambda))/sqrt(2))/sqrt( abs(alpha)^2+abs(lambda)^2+abs(ro)^2 );
    
    Orientation=(phase(ro)-phase(lambda))/4 ;
    R=0.5*[cos(2*Orientation)+1 -sin(2*Orientation) -sin(2*Orientation) 1-cos(2*Orientation);
    sin(2*Orientation) cos(2*Orientation)+1 cos(2*Orientation)-1 -sin(2*Orientation);
    sin(2*Orientation) cos(2*Orientation)-1 cos(2*Orientation)+1 -sin(2*Orientation);
    1-cos(2*Orientation) sin(2*Orientation) sin(2*Orientation) cos(2*Orientation)+1];
    
    Ori_max(kk)=Orientation*(180/pi);%atand(sin(Orientation)/cos(Orientation));
    Sym_max=R*(alpha*(1/sqrt(2))*[1;0;0;1]+(exp(1i*(phase(ro)+phase(lambda))/2)/sqrt(2))*(abs(lambda)+abs(ro))*(1/sqrt(2))*[1;0;0;-1]);

    
    % Parameterization
    L1=a*cos_taw*(1/norm(Sym_max))*alpha;
    L2=a*cos_taw*(1/norm(Sym_max))*(exp(1i*(phase(ro)+phase(lambda))/2)/sqrt(2))*(abs(lambda)+abs(ro));
    L3=a*sin_taw*(-1i)*exp(1i*(phase(ro)+phase(lambda))/2);

    Alphas(kk) = atan2d(abs(L2),abs(L1));
    Phis(kk) = angle(L2*conj(L1))*180/pi;
    Taw(kk) = asind(  abs(L3)./sqrt(abs(L3)^2+abs(L2)^2+abs(L1)^2));
    M(kk)= sqrt(abs(L3)^2+abs(L2)^2+abs(L1)^2);

    %______________________________________________________________________
    % disp 
%         parfor_progress;
    if rem(kk, dispN)==0
       kk    
    end
    %______________________________________________________________________
end

Alphas = reshape(Alphas,Lrng,Lazi)';
Phis = reshape(Phis,Lrng,Lazi)';
Taw = reshape(Taw,Lrng,Lazi)';
Ori_max = reshape(Ori_max,Lrng,Lazi)';
M = reshape(M,Lrng,Lazi)';
