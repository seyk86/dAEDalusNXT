%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Lift1,Pitching_Moment1,c_lift,c_moment,Q_theodorsen]=Theodorsen_Crepaldi2(rho,a,b,c,Uinf,hh,ah,k)

% function to calculate aerodynamic coefficients based on Theodorsen unsteady method 

% 29.07.2016 - Error: dimensionless coefficients changing with the velocity

%% Input parameters
% rho: air density
% a: pitch axis, and reference axis for pitching moment coordinate
% b: semi-chord
% c: hinge axis coordinate
% Uinf: freestream velocity
% hh: heave amplitude
% ah: pitch amplitude
% k: reduced frequency

%% Definitions
% a=-.4; c=0.6; b=1; % a is the COORDINATE (origin on the semichord) of the pitch axis for pitching moment (+TE,-LE)
chord=2*b;

sqr=sqrt(1-c*c);
arc=acos(c);
heave_amplitude=hh;
pitch_amplitude=ah;
time=0:0.001:10;


%% T functions
t1=-(2+c*c)/3*sqr+c*arc;
t3=-(1-c*c)/8*(5.*c*c+4)+.25*c*(7+2*c*c)*sqr*arc-(1./8.+c*c)*arc*arc;
t4=c*sqr-arc;
t5=-(1-c*c)-arc*arc+2*c*sqr*arc;
t7=c*(7+2*c*c)/8*sqr-(1/8+c*c)*arc;
t8=-1/3*(1+2*c*c)*sqr+c*arc;
t9=.5*(sqr*(1-c*c)/3+a*t4);
t10=sqr+arc;
t11=(2-c)*sqr+(1-2*c)*arc;
t12=(2+c)*sqr-(1+2*c)*arc;
t13=-.5*(t7+(c-a)*t1);

%% Useful T functions
t15=t4+t10;
t16=t1-t8-(c-a)*t4+0.5*t11;
t17=-2*t9-t1+(a-0.5)*t4;
t18=t5-t4*t10;
t19=-0.5*t4*t11;

%% Preallocation
Lift1=zeros(2,size(k,2));
Pitching_Moment1=Lift1;
Lift2=Lift1;
Pitching_Moment2=Lift1;
Lift2_1=Lift1;
Lift2_2=Lift1;
Lift3=Lift1;
Pitching_Moment3=Lift1;
Lift4=Lift1;

%% Loop

for Movement=1:2
    
    
    %% Running Heave and Pitch movement
    if Movement==1 %Heave
        h_0=heave_amplitude;
        alpha_0=0;
        fprintf(['Theodorsen for HEAVE oscillation:                 ' '\n']);
    elseif Movement==2 %Pitch
        h_0=0;
        alpha_0=pitch_amplitude;
        fprintf(['Theodorsen for PITCH oscillation:                 ' '\n']);
    else %Heave and Pitch
        h_0=heave_amplitude;
        alpha_0=pitch_amplitude;
        fprintf(['Theodorsen for HEAVE and PITCH oscillation:       ' '\n']);
    end
    
    
    
    %% Theodorsen function for each reduced frequency
    for j=1:size(k,2)
        
        %% Input parameters
        
        k_j=k(j);
        omega=Uinf*k_j/b;
        %         angle_max(j,1)=atand(2*heave_amplitude*k/chord);
        %         angle_max(j,2)=2*pi*angle_max(j,1)*pi/180;
        fprintf('\b\b\b\b\b %02d',round(j/size(k,2)*100));
        fprintf(['%%' '\n']);
        
        
        
        %% Time-dependent parameters(h, alpha) and cosine curve fitting
        
        warning('off','all')
        
        h=h_0*sin(omega*time);
        h_dot=omega*h_0*cos(omega*time);
        h_2dot=-omega^2*h;
        
        x_h=Lsfit_osc(omega,time,h);
        x_h_dot=Lsfit_osc(omega,time,h_dot);
        x_h_2dot=Lsfit_osc(omega,time,h_2dot);
        
        alpha=alpha_0*sin(omega*time);
        alpha_dot=omega*alpha_0*cos(omega*time);
        alpha_2dot=-omega^2*alpha;
        
        x_alpha=Lsfit_osc(omega,time,alpha);
        x_alpha_dot=Lsfit_osc(omega,time,alpha_dot);
        x_alpha_2dot=Lsfit_osc(omega,time,alpha_2dot);
        
        beta=zeros(size(alpha));
        beta_dot=zeros(size(alpha_dot));
        beta_2dot=zeros(size(alpha_2dot));
        
        x_beta=Lsfit_osc(omega,time,beta);
        x_beta_dot=Lsfit_osc(omega,time,beta_dot);
        x_beta_2dot=Lsfit_osc(omega,time,beta_2dot);
        
        if k_j>0
            input_vec=[x_h(1)*exp(1i*(x_h(2)+pi/2)) ;
                x_alpha(1)*exp(1i*(x_alpha(2)+pi/2)) ;
                x_beta(1)*exp(1i*(x_beta(2)+pi/2))]; % h, alpha, beta inputs
            
            input_vec_dot=[x_h_dot(1)*exp(1i*(x_h_dot(2)+pi/2)) ;
                x_alpha_dot(1)*exp(1i*(x_alpha_dot(2)+pi/2)) ;
                x_beta_dot(1)*exp(1i*(x_beta_dot(2)+pi/2))]; % h, alpha, beta inputs (first order derivative)
            
            input_vec_2dot=[x_h_2dot(1)*exp(1i*(x_h_2dot(2)+pi/2)) ;
                x_alpha_2dot(1)*exp(1i*(x_alpha_2dot(2)+pi/2)) ;
                x_beta_2dot(1)*exp(1i*(x_beta_2dot(2)+pi/2))]; % h, alpha, beta inputs (second order derivative)
        else
            input_vec=zeros(3,1);
            input_vec_dot=input_vec;
            input_vec_2dot=input_vec;
            
            input_vec(1,1)=h_0;
            input_vec(2,1)=alpha_0;
        end
        
        
        
        %% Parameters for the Theodorsen formula (kinematic parameters in complex domain and theirs derivatives)
        h=input_vec(1);
        h1=input_vec_dot(1);
        h2=input_vec_2dot(1);
        alpha=input_vec(2);
        alpha1=input_vec_dot(2);
        alpha2=input_vec_2dot(2);
        beta=input_vec(3);
        beta1=input_vec_dot(3);
        beta2=input_vec_2dot(3);
        
        
        
        %% Calculation
        
        %Lift deficiency factor (choose function)
        
        %         [ck_imag,ck_real,ck_magnit,ck_phase]=Theodorsen_parameters(k_j);
        [ck_imag,ck_real,ck_magnit,ck_phase]=Theodorsen_parameters2(k_j);
        C=ck_real+1i*ck_imag;
        
        %Common used formula (from Bisplinghoff book)
        Lift1(Movement,j)=pi*rho*b^2*(h2+Uinf*alpha1-b*a*alpha2)+... % non-circulatory
            2*pi*rho*Uinf*b*C*(h1+Uinf*alpha+b*(0.5-a)*alpha1); % circulatory
        Pitching_Moment1(Movement,j)=pi*rho*b^2*(b*a*h2-...
            Uinf*b*(0.5-a)*alpha1-b^2*((1/8)+a^2)*alpha2)+...
            2*pi*rho*Uinf*b^2*(a+0.5)*(h1+Uinf*alpha+b*(0.5-a)*alpha1)*C;
        
        %Katz
        Lift2(Movement,j)=pi*rho*Uinf*chord*C*(Uinf*alpha-h1+(.75-(a/chord))*chord*alpha1)+...
            (pi*rho*chord^2/4)*((Uinf*alpha1-h2)+chord*(.5-(a/chord))*alpha2);
        Pitching_Moment2(Movement,j)=chord*h2*((a/chord)-.5)+...
            Uinf*chord*alpha1*(.75-(a/chord))+...
            chord^2/4*((9/8)+4*(a^2/chord^2)-4*(a/chord))*alpha2-...
            ((4*a/chord)-1)*Uinf*C*(-h1+Uinf*alpha+chord*alpha1*(.75-(a/chord)));
        Pitching_Moment2(Movement,j)=-(pi*rho*chord^2/4)*Pitching_Moment2(Movement,j);
        Lift2_1(Movement,j)=pi*rho*Uinf*chord*C*(Uinf*alpha-h1+(.75-(a/chord))*chord*alpha1);
        Lift2_2(Movement,j)= (pi*rho*chord^2/4)*((Uinf*alpha1-h2)+chord*(.5-(a/chord))*alpha2);
        
        %Theodorsen paper
        Lift3(Movement,j)=-rho*b^2*(Uinf*pi*alpha1+pi*h2-pi*b*a*alpha2-Uinf*t4*beta1-t1*b*beta2)-...
            2*pi*rho*Uinf*b*C*(Uinf*alpha+h1+b*(0.5-a)*alpha1+(1/pi)*t10*Uinf*beta+b*(1/(2*pi))*t11*beta1);
        Pitching_Moment3(Movement,j)=-rho*b^2*(pi*(0.5-a)*Uinf*b*alpha1+pi*b^2*((1/8)+a^2)*alpha2+...
            (t4+t10)*Uinf^2*beta+...
            (t1-t8-(c-a)*t4+0.5*t11)*Uinf*b*beta1-...
            (t7+(c-a)*t1)*b^2*beta2-a*pi*b*h2)+...
            2*rho*Uinf*b^2*pi*(a+0.5)*C*(Uinf*alpha+h1+b*(0.5-a)*alpha1+...
            (1/pi)*t10*Uinf*beta+b*(1/2/pi)*t11*beta1);
        
        %Brunton-theodorsen
        Lift4(Movement,j)=pi*(h2+h1-a*alpha2)+...
            2*pi*(alpha+h1+alpha1*(0.5-a))*C;
        
    end
    
end




%% Normalizations

dyn_press=0.5*rho*Uinf^2;

%Choose method
Lift=-Lift1;
Pitching_Moment=Pitching_Moment1;

%Dimensionless parameters
c_lift=-Lift/(dyn_press*chord);
c_moment=Pitching_Moment/(dyn_press*chord^2);

Lift_norm(1,:)=Lift(1,:)*(1/heave_amplitude);
Lift_norm(2,:)=Lift(2,:)*(1/pitch_amplitude);
Pitching_Moment_norm(1,:)=Pitching_Moment(1,:)*(1/heave_amplitude);
Pitching_Moment_norm(2,:)=Pitching_Moment(2,:)*(1/pitch_amplitude);




%% Time Plot
% fprintf(['Plotting' '\n'])

% for Movement=1:1
% %     figure;
%     figure(1)
%     for jjj=1:size(k,2)
%         plot(time,abs(c_lift(Movement,jjj))*sin((Uinf*k(jjj)/b)*time+angle(c_lift(Movement,jjj))));
% %         plot(time,abs(Lift2_1(Movement,jjj))*sin((Uinf*k(jjj)/b)*time-ck_phase)+...
% %             abs(Lift2_2(Movement,jjj))*sin((Uinf*k(jjj)/b)*time));
%         xlabel('Time')
%         ylabel('Lift coefficient')
%         hold on
%         grid on
%         legend_aux = cellstr(num2str(k', 'k=%.3g'));
%         legend(legend_aux)
%     end
% end
%
% for Movement=1:1
% %     figure;
%     figure(2)
%     for jjj=1:size(k,2)
%         plot(time,abs(c_moment(Movement,jjj))*sin((Uinf*k(jjj)/b)*time+angle(c_moment(Movement,jjj))));
%         xlabel('Time')
%         ylabel('Pitching moment coefficient')
%         hold on
%         grid on
%         legend_aux = cellstr(num2str(k', 'k=%.3g'));
%         legend(legend_aux)
%     end
% end



%% Klaus comparison

% Lift=-Lift1/(dyn_press);
% Pitching_Moment=Pitching_Moment1/(dyn_press);
%
% Lift(1,:)=Lift(1,:)*(1/heave_amplitude);
% Pitching_Moment(1,:)=Pitching_Moment(1,:)*(1/heave_amplitude);
% Lift(2,:)=Lift(2,:)*(1/pitch_amplitude);
% Pitching_Moment(2,:)=Pitching_Moment(2,:)*(1/(pitch_amplitude));
%
% figure;
%
% subplot(2,2,1)
% plot(k,imag(Lift(1,:)),'b--*')
% hold on
% plot(k,real(Lift(1,:)),'r--*')
% grid on
% legend('Imag','Real')
% title('Lift Heave')
% xlabel('Reduced frequency')
%
% subplot(2,2,2)
% plot(k,imag(Lift(2,:)),'b--*')
% hold on
% plot(k,real(Lift(2,:)),'r--*')
% grid on
% legend('Imag','Real')
% title('Lift Pitch')
% xlabel('Reduced frequency')
%
% subplot(2,2,3)
% plot(k,imag(Pitching_Moment(1,:)),'b--*')
% hold on
% plot(k,real(Pitching_Moment(1,:)),'r--*')
% grid on
% legend('Imag','Real')
% title('Pitching Moment Heave')
% xlabel('Reduced frequency')
%
% subplot(2,2,4)
% plot(k,imag(Pitching_Moment(2,:)),'b--*')
% hold on
% plot(k,real(Pitching_Moment(2,:)),'r--*')
% grid on
% legend('Imag','Real')
% title('Pitching Moment Pitch')
% xlabel('Reduced frequency')



%% Save parameters as dAEDalus
% fprintf(['Saving' '\n'])


% Q_theodorsen(:,1,kk)=Lift_norm(:,1)/dyn_press;
% Q_theodorsen(:,2,kk)=Pitching_Moment_norm(:,1)/dyn_press;

Q_theodorsen(:,1,size(k,2))=Lift_norm(:,1)/dyn_press;
Q_theodorsen(:,2,size(k,2))=Pitching_Moment_norm(:,1)/dyn_press;

% cd('C:\Users\Studer\Desktop\Unsteady_Replicated_Values_Klaus')
% save('Q_implementation_Theodorsen2.mat','Q_theodorsen')
% cd('C:\Users\Studer\Desktop\Unsteady_Replicated_Values_Klaus')
% Timee=toc;
% fprintf(['Time: ' num2str(Timee) '\n'])

end
