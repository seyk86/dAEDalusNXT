%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [A0,A1,A2,E,D]=minimum_state_approximation(Q,k,gamma,nits) 
%% rogers approximation
offset=0;

nQ=length(Q(:,1,1));
nL=length(gamma);
A_base=zeros(2*length(k)*nQ,3+nL);

Wr=ones(nQ,nQ);
Wi=ones(nQ,nQ);

%data normalization
 epsi=1;
 for kk=1:length(k)
     for i=1:nQ
         for j=1:nQ
             max_r=max(max(real(Q(i,j,:))),epsi);
             max_i=max(max(imag(Q(i,j,:))),epsi);
             Wr(i,j,kk)=epsi/(max_r*k(kk));
             Wi(i,j,kk)=epsi/(max_i*k(kk));
         end
     end
 end
% 
% for i=1:nQ
%     for j=1:nQ
%         max_r=max(max(real(Q(i,j,:))),epsi);
%         max_i=max(max(imag(Q(i,j,:))),epsi);
%         if max_r<=epsi
%             Wr(i,j)=1;
%         else
%             Wr(i,j)=1/max_r;
%         end
%         if max_i<=epsi
%             Wi(i,j)=1;
%         else
%             Wi(i,j)=1/max_r;
%         end
%     end
% end

% Wr(7:8,7:8)=2.3; % wouldn't be correct for the nonlinear optimization
% Wi(7:8,7:8)=2.3;

% Wr(7:end,7:end)=120;
% Wi(7:end,7:end)=120;
% Wr(7:8,7:8)=280;
% Wi(7:8,7:8)=280;

E=zeros(nL,nQ);
D=ones(nQ,nL);
D(1:min(nL,nQ),1:min(nL,nQ))=eye(min(nL,nQ),min(nL,nQ));

for interations=1:1:nits
    redi=0;
    for col_loop=1:1:nQ
        offset=0;
        for ii=1:1:nQ
            k_idx=1;
            for i=1:2:2*length(k)    
                b(i+offset)=real(Q(ii,col_loop,k_idx))*Wr(ii,col_loop,ceil(i/2));
                b(i+1+offset)=imag(Q(ii,col_loop,k_idx))*Wi(ii,col_loop,ceil(i/2));
                A_base(i+offset,1+3*(ii-1):3*ii)=[1 0 -k(k_idx)^2]*Wr(ii,col_loop,ceil(i/2));
                A_base(i+1+offset,1+3*(ii-1):3*ii)=[0  k(k_idx) 0]*Wi(ii,col_loop,ceil(i/2));
                for j=1:1:nL
                    A_base(i+offset,3*nQ+j)=D(ii,j)*k(k_idx)^2/(k(k_idx)^2+gamma(j)^2)*Wr(ii,col_loop,ceil(i/2));
                    A_base(i+offset+1,3*nQ+j)=D(ii,j)*k(k_idx)*gamma(j)/(k(k_idx)^2+gamma(j)^2)*Wi(ii,col_loop,ceil(i/2));
                end
                k_idx=k_idx+1;
            end
            offset=offset+length(k)*2;
        end
        %% Solution for Rogers Approximation#
        [sol,nn,resid]=lsqlin(A_base,b');
                redi=redi+sum(abs(resid));
        A0(:,col_loop)=sol(1:3:3*nQ);
        A1(:,col_loop)=sol(2:3:3*nQ);
        A2(:,col_loop)=sol(3:3:3*nQ);
        E(1:nL,col_loop)=sol(3*nQ+1:end);
    end
     redi   
    % D computation
    A_base=A_base*0;
    b=b*0;
    for row_loop=1:1:nQ
        offset=0;
        for ii=1:1:nQ
            k_idx=1;
            for i=1:2:2*length(k)
                b(i+offset)=real(Q(row_loop,ii,k_idx))*Wr(row_loop,ii,ceil(i/2));
                b(i+1+offset)=imag(Q(row_loop,ii,k_idx))*Wi(row_loop,ii,ceil(i/2));
                A_base(i+offset,1+3*(ii-1):3*ii)=[1 0 -k(k_idx)^2]*Wr(row_loop,ii,ceil(i/2));
                A_base(i+offset+1,1+3*(ii-1):3*ii)=[0    k(k_idx) 0]*Wi(row_loop,ii,ceil(i/2));
                for j=1:1:nL
                    A_base(i+offset,3*nQ+j)=E(j,ii)*k(k_idx)^2/(k(k_idx)^2+gamma(j)^2)*Wr(row_loop,ii,ceil(i/2));
                    A_base(i+offset+1,3*nQ+j)=E(j,ii)*k(k_idx)*gamma(j)/(k(k_idx)^2+gamma(j)^2)*Wi(row_loop,ii,ceil(i/2));
                end
                k_idx=k_idx+1;
            end
            offset=offset+length(k)*2;
        end
        %% Solution for Rogers Approximation
        sol=lsqlin(A_base,b');
        A0(row_loop,:)=sol(1:3:3*nQ);
        A1(row_loop,:)=sol(2:3:3*nQ);
        A2(row_loop,:)=sol(3:3:3*nQ);
        D(row_loop,1:nL)=sol(3*nQ+1:end);
    end
    A_base=A_base*0;
    b=b*0;  

end

    
    nQ=size(Q,1);
    debug=1;
    if debug==1
        for i=1:nQ
            for j=1:nQ
                kk=1;
                for kc=min(k):0.01:max(k)
                    Qtilde(kk,1)=A0(i,j)-kc^2*A2(i,j);
                    Qtilde(kk,2)=kc*A1(i,j);
                    for lag=1:nL
                        Qtilde(kk,1)=Qtilde(kk,1)+E(lag,j)*D(i,lag)*kc^2/(kc^2+gamma(lag)^2);
                    end
                    
                    for lag=1:nL
                        Qtilde(kk,2)=Qtilde(kk,2)+E(lag,j)*D(i,lag)*kc*gamma(lag)/(kc^2+gamma(lag)^2);
                    end
                    kk=kk+1;
                end
                %figure
                %plot(real(squeeze(Q(i,j,:))),imag(squeeze(Q(i,j,:))),'-rx');
                %hold on
                %plot(Qtilde(:,1),Qtilde(:,2),'-b');
                %legend('Data','Minimum State Approximation (N=6)')
            end
        end
    end
