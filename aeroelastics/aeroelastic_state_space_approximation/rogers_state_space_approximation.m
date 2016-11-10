%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [A0, A1, A2, Arest] = rogers_state_space_approximation(Q,k,gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nL=length(gamma);
nQ=size(Q,1);

dim_F=2+nL;
dim_G=1+nL;
dim_A=dim_F+dim_G;
for i=1:nQ
    for j=1:nQ
        k_idx=1;
        D=ones(1,nL+3);
        E=zeros(nL+1,1);
        W=eye(nL+3,nL+3);
        A_base=zeros(2*length(k),3+nL);
        for kk=1:2:2*length(k)
            b(kk)=real(Q(i,j,k_idx));
            b(kk+1)=imag(Q(i,j,k_idx));
            A_base(kk,1:3)=[1 0 -k(k_idx)^2];
            A_base(kk+1,1:3)=[0  k(k_idx) 0];
            for lag=1:1:nL
                A_base(kk,3+lag)=k(k_idx)^2/(k(k_idx)^2+gamma(lag)^2);
                A_base(kk+1,3+lag)=k(k_idx)*gamma(lag)/(k(k_idx)^2+gamma(lag)^2);
            end
            k_idx=k_idx+1;
        end
        
        %% Solution for Rogers Approximation
        %E=A_base\b';
        E=lsqlin(A_base,b);
        E_full(i,j)=E(1);
        E_full(i+nQ,j)=E(2);
        E_full(i+2*nQ,j)=E(3);
        for lag=1:nL
            E_full(i+(2+lag)*nQ,j)=E(3+lag);
        end
    end
end

%E_full=X_roger;

A0=E_full(1:nQ,1:nQ);
A1=E_full(1+nQ:2*nQ,1:nQ);
A2=E_full(1+2*nQ:3*nQ,1:nQ);

%rearrage AN
for lag=1:nL
    Arest(1:nQ,1+(lag-1)*nQ:lag*nQ)=E_full(1+(2+lag)*nQ:(3+lag)*nQ,1:nQ);
end

debug=0;

if debug==1
    for i=1:nQ
        j=5;
        %for j=1:nQ
            kk=1;
            for kc=0:0.01:max(k)
                Qtilde(kk,1)=0;
                Qtilde(kk,2)=0;

                Qtilde(kk,1)=A0(i,j)-kc^2*A2(i,j);
                Qtilde(kk,2)=kc*A1(i,j);
                for lag=1:nL
                  %  Qtilde(kk,1)=Qtilde(kk,1)+kc^2*E_full(i+(2+lag)*nQ,j)/(kc^2+gamma(lag)^2);
                    Qtilde(kk,1)=Qtilde(kk,1)+kc^2*Arest(i,j+(lag-1)*nQ)/(kc^2+gamma(lag)^2);
                end
                
                for lag=1:nL
                   % Qtilde(kk,2)=Qtilde(kk,2)+kc*E_full(i+(2+lag)*nQ,j)*gamma(lag)/(kc^2+gamma(lag)^2);
                    Qtilde(kk,2)=Qtilde(kk,2)+kc*Arest(i,j+(lag-1)*nQ)*gamma(lag)/(kc^2+gamma(lag)^2);
                end
                kk=kk+1;
            end
            figure
            plot(real(squeeze(Q(i,j,:))),imag(squeeze(Q(i,j,:))),'-rx');
            hold on
            plot(Qtilde(:,1),Qtilde(:,2),'-bx');
            legend('Data','Rogers Approximation (N=6)')
       % end
    end
end

end

