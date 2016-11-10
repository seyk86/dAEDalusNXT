%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   This file is part of dAEDalus structures
%                   Copyright (C) 2011, Klaus Seywald
%     Author:   	Klaus Seywald
%                   klaus.seywald@mytum.de
%                   seywald@kth.se
% lin_elK_6dof:     class method of class_beamelement
%                   calculates 6DOF linear stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj=lin_elK_6dof(obj)

    E=obj.E;
    G=obj.G;
    A=obj.A;
    L=obj.le;
    Ix=obj.Ix;
    Iz=obj.Iz;
    Ip=obj.Ip;
    Izx=obj.Izx;
    EIz=E*Iz;
    EIzx=E*Izx;
    EIx=E*Ix;
    GJ=G*Ip;
    EA=E*A;
%    Calculate Element Stiffness Matrices in Local Coordinates
    obj.elK=[12*EIz/L^3  ,0     ,0  ,0  ,0  ,-6*EIz/L^2  ,-12*EIz/L^3  ,0      ,0      ,0  ,0  ,-6*EIz/L^2;...
             0           ,EA/L  ,0  ,0  ,0  ,0           ,0            ,-EA/L  ,0  ,0  ,0,0;...
        0,0,12*EIx/L^3,6*EIx/L^2,0,0,0,0,-12*EIx/L^3,6*EIx/L^2,0,0;...
        0,0,6*EIx/L^2,4*EIx/L,0,0,0,0,-6*EIx/L^2,2*EIx/L,0,0;...
        0,0,0,0,GJ/L,0,0,0,0,0,-GJ/L,0;...
        -6*EIz/L^2,0,0,0,0,4*EIz/L,6*EIz/L^2,0,0,0,0,2*EIz/L;...
        -12*EIz/L^3,0,0,0,0,6*EIz/L^2,12*EIz/L^3,0,0,0,0,6*EIz/L^2;...
        0,-EA/L,0,0,0,0,0,EA/L,0,0,0,0;...
        0,0,-12*EIx/L^3,-6*EIx/L^2,0,0,0,0,12*EIx/L^3,-6*EIx/L^2,0,0;...
        0,0,6*EIx/L^2,2*EIx/L,0,0,0,0,-6*EIx/L^2,4*EIx/L,0,0;...
        0,0,0,0,-GJ/L,0,0,0,0,0,GJ/L,0;...
        -6*EIz/L^2,0,0,0,0,2*EIz/L,6*EIz/L^2,0,0,0,0,4*EIz/L;...
        ];
    %% New local stiffness matrix definition
    %BEFORE THE CALCULATION, EIzx, AT, LT and FT needs to be defined.
    %AT=coupled axial-torsion stiffness
    %FT=coupled flap-torsion stiffness
    %LT=coupled lag-torsion stiffness
    obj.elK=zeros(12);
    
    LT=0;
    AT=0;
    FT=0;
    
    
%     
     %% ORIGINAL
%     %OBS! The output order is
%     %uy1,ux1,uz1,fi1,tetax1,tetaz1,uy2,ux2,uz2,fi2,tetax2,tetaz2. So the
%     %lines and collumns needs to be rearrenged.
%     %First block (1,1) to (6,6)
%     obj.elK(1,1)=EA/L  ;obj.elK(1,2)=0           ;obj.elK(1,3)=0                 ;obj.elK(1,4)=AT/L          ;obj.elK(1,5)=0              ;obj.elK(1,6)=0            ; %OK
%                         obj.elK(2,2)=12*EIz/L^3  ;obj.elK(2,3)=12*EIzx/L^3       ;obj.elK(2,4)=0             ;obj.elK(2,5)=-6*EIzx/L^2    ;obj.elK(2,6)=-6*EIz/L^2   ; %OK
%                                                   obj.elK(3,3)=12*EIx/L^3        ;obj.elK(3,4)=0             ;obj.elK(3,5)=6*EIx/L^2      ;obj.elK(3,6)=6*EIzx/L^2   ; %OK
%                                                                                   obj.elK(4,4)=GJ/L          ;obj.elK(4,5)=FT/L           ;obj.elK(4,6)=-LT/L        ; %OK
%                                                                                                               obj.elK(5,5)=4*EIx/L        ;obj.elK(5,6)=-4*EIzx/L    ; %OK
%                                                                                                                                            obj.elK(6,6)=4*EIz/L      ; %OK
%                                                                                                                                            
%     %Second block (1,7) to (6,12)
%     obj.elK(1,7)=-obj.elK(1,1);  obj.elK(1,8)=-obj.elK(1,2);  obj.elK(1,9)=-obj.elK(1,3);  obj.elK(1,10)=-obj.elK(1,4);  obj.elK(1,11)=-obj.elK(1,5);  obj.elK(1,12)=-obj.elK(1,6);%OK
%     obj.elK(2,7)=-obj.elK(1,2);  obj.elK(2,8)=-obj.elK(2,2);  obj.elK(2,9)=-obj.elK(2,3);  obj.elK(2,10)=-obj.elK(2,4);  obj.elK(2,11)=+obj.elK(2,5);  obj.elK(2,12)=+obj.elK(2,6);%OK
%     obj.elK(3,7)=-obj.elK(1,3);  obj.elK(3,8)=-obj.elK(2,3);  obj.elK(3,9)=-obj.elK(3,3);  obj.elK(3,10)=-obj.elK(3,4);  obj.elK(3,11)=+obj.elK(3,5);  obj.elK(3,12)=+obj.elK(3,6);%OK
%     obj.elK(4,7)=-obj.elK(1,4);  obj.elK(4,8)=-obj.elK(2,4);  obj.elK(4,9)=-obj.elK(3,4);  obj.elK(4,10)=-obj.elK(4,4);  obj.elK(4,11)=-FT/L        ;  obj.elK(4,12)=LT/L         ;%OK
%     obj.elK(5,7)=-obj.elK(1,5);  obj.elK(5,8)=-obj.elK(2,5);  obj.elK(5,9)=-obj.elK(3,5);  obj.elK(5,10)=-FT/L        ;  obj.elK(5,11)=2*EIx/L      ;  obj.elK(5,12)=-2*EIzx/L    ;%OK
%     obj.elK(6,7)=-obj.elK(1,6);  obj.elK(6,8)=-obj.elK(2,6);  obj.elK(6,9)=-obj.elK(3,6);  obj.elK(6,10)=LT/L         ;  obj.elK(6,11)=-2*EIzx/L    ;  obj.elK(6,12)=2*EIz/L      ;%OK
%    
%     
%     %Third block (7,7) to (12,12)
%     obj.elK(7,7)=obj.elK(1,1);   obj.elK(7,8)=0           ;   obj.elK(7,9)=0           ; obj.elK(7,10)= AT/L        ;  obj.elK(7,11)=0            ;  obj.elK(7,12)=0            ;%OK
%                                  obj.elK(8,8)=obj.elK(2,2);   obj.elK(8,9)=obj.elK(2,3); obj.elK(8,10)=obj.elK(2,4) ;  obj.elK(8,11)=-obj.elK(2,5);  obj.elK(8,12)=-obj.elK(2,6);%OK
%                                                               obj.elK(9,9)=obj.elK(3,3); obj.elK(9,10)=obj.elK(3,4) ;  obj.elK(9,11)=-obj.elK(3,5);  obj.elK(9,12)=-obj.elK(3,6);%OK
%                                                                                          obj.elK(10,10)=obj.elK(4,4);  obj.elK(10,11)=obj.elK(4,5);  obj.elK(10,12)=obj.elK(4,6);%OK
%     
%                       
%                                                                                          
%                                        
%                                                                                          
%                                                                                          
                                                                                         
    %% SIGNS MODIFIED                                                                                     
    obj.elK(1,1)=EA/L  ;obj.elK(1,2)=0           ;obj.elK(1,3)=0                 ;obj.elK(1,4)=AT/L          ;obj.elK(1,5)=0              ;obj.elK(1,6)=0            ; %OK
                        obj.elK(2,2)=12*EIz/L^3  ;obj.elK(2,3)=12*EIzx/L^3       ;obj.elK(2,4)=0             ;obj.elK(2,5)=-6*EIzx/L^2    ;obj.elK(2,6)=-6*EIz/L^2   ; %OK   %obj.elK(2,6)=-6*EIz/L^2 (same as original matrix)
                                                  obj.elK(3,3)=12*EIx/L^3        ;obj.elK(3,4)=0             ;obj.elK(3,5)=6*EIx/L^2     ;obj.elK(3,6)=6*EIzx/L^2   ; %OK    %obj.elK(3,5)=6*EIx/L^2 (same as original)
                                                                                  obj.elK(4,4)=GJ/L          ;obj.elK(4,5)=FT/L           ;obj.elK(4,6)=-LT/L        ; %OK
                                                                                                              obj.elK(5,5)=4*EIx/L        ;obj.elK(5,6)=-4*EIzx/L    ; %OK
                                                                                                                                           obj.elK(6,6)=4*EIz/L      ; %OK
                                                                                                                                           
    %Second block (1,7) to (6,12)
    obj.elK(1,7)=-obj.elK(1,1);  obj.elK(1,8)=-obj.elK(1,2);  obj.elK(1,9)=-obj.elK(1,3);  obj.elK(1,10)=-obj.elK(1,4);  obj.elK(1,11)=-obj.elK(1,5);  obj.elK(1,12)=-obj.elK(1,6);%OK
    obj.elK(2,7)=-obj.elK(1,2);  obj.elK(2,8)=-obj.elK(2,2);  obj.elK(2,9)=-obj.elK(2,3);  obj.elK(2,10)=-obj.elK(2,4);  obj.elK(2,11)=+obj.elK(2,5);  obj.elK(2,12)=+obj.elK(2,6);%OK
    obj.elK(3,7)=-obj.elK(1,3);  obj.elK(3,8)=-obj.elK(2,3);  obj.elK(3,9)=-obj.elK(3,3);  obj.elK(3,10)=-obj.elK(3,4);  obj.elK(3,11)=+obj.elK(3,5);  obj.elK(3,12)=+obj.elK(3,6);%OK
    obj.elK(4,7)=-obj.elK(1,4);  obj.elK(4,8)=-obj.elK(2,4);  obj.elK(4,9)=-obj.elK(3,4);  obj.elK(4,10)=-obj.elK(4,4);  obj.elK(4,11)=-FT/L        ;  obj.elK(4,12)=LT/L         ;%OK
    obj.elK(5,7)=-obj.elK(1,5);  obj.elK(5,8)=-obj.elK(2,5);  obj.elK(5,9)=-obj.elK(3,5);  obj.elK(5,10)=-FT/L        ;  obj.elK(5,11)=2*EIx/L      ;  obj.elK(5,12)=-2*EIzx/L    ;%OK
    obj.elK(6,7)=-obj.elK(1,6);  obj.elK(6,8)=-obj.elK(2,6);  obj.elK(6,9)=-obj.elK(3,6);  obj.elK(6,10)=LT/L         ;  obj.elK(6,11)=-2*EIzx/L    ;  obj.elK(6,12)=2*EIz/L      ;%OK
   
    
    %Third block (7,7) to (12,12)
    obj.elK(7,7)=obj.elK(1,1);   obj.elK(7,8)=0           ;   obj.elK(7,9)=0           ; obj.elK(7,10)= AT/L        ;  obj.elK(7,11)=0            ;  obj.elK(7,12)=0            ;%OK
                                 obj.elK(8,8)=obj.elK(2,2);   obj.elK(8,9)=obj.elK(2,3); obj.elK(8,10)=obj.elK(2,4) ;  obj.elK(8,11)=-obj.elK(2,5);  obj.elK(8,12)=-obj.elK(2,6);%OK
                                                              obj.elK(9,9)=obj.elK(3,3); obj.elK(9,10)=obj.elK(3,4) ;  obj.elK(9,11)=-obj.elK(3,5);  obj.elK(9,12)=-obj.elK(3,6);%OK
                                                                                         obj.elK(10,10)=obj.elK(4,4);  obj.elK(10,11)=obj.elK(4,5);  obj.elK(10,12)=obj.elK(4,6);%OK
                                                                                                                       obj.elK(11,11)=obj.elK(5,5);  obj.elK(11,12)=obj.elK(5,6);%OK
                                                                                                                                                     obj.elK(12,12)=obj.elK(6,6);%OK                                                                                                                   obj.elK(11,11)=obj.elK(5,5);  obj.elK(11,12)=obj.elK(5,6);%OK
                                                                                                                                                     obj.elK(12,12)=obj.elK(6,6);%OK
    %Defining the simmetric matrix
    obj.elK=obj.elK+tril(obj.elK.',-1);
    
    %Now that we have the simmetric matrix, we need to relocate their
    %collumns and lines to adjust the inputs and outputs order
    % 1<->2
    % 4<->5
    % 7<->8
    % 10<->11
    %To do that, we need to save the lines and collumns that needs to be
    %changed.
    Savematrix=obj.elK;
    
    %Change of collumn
    for j=1:1:12
        if j==1
            obj.elK(:,j+1)=Savematrix(:,j);
        elseif j==4
            obj.elK(:,j+1)=Savematrix(:,j);
        elseif j==7
            obj.elK(:,j+1)=Savematrix(:,j);
        elseif j==10
            obj.elK(:,j+1)=Savematrix(:,j);
        end
        if j==2
            obj.elK(:,j-1)=Savematrix(:,j);
        elseif j==5
            obj.elK(:,j-1)=Savematrix(:,j);
        elseif j==8
            obj.elK(:,j-1)=Savematrix(:,j);
        elseif j==11
            obj.elK(:,j-1)=Savematrix(:,j);
        end
    end
    Savematrix1=obj.elK;
    %Change of lines
    for j=1:1:12
        if j==1
            obj.elK(j+1,:)=Savematrix1(j,:);
        elseif j==4
            obj.elK(j+1,:)=Savematrix1(j,:);
        elseif j==7
            obj.elK(j+1,:)=Savematrix1(j,:);
        elseif j==10
            obj.elK(j+1,:)=Savematrix1(j,:);
        end
        if j==2
            obj.elK(j-1,:)=Savematrix1(j,:);
        elseif j==5
            obj.elK(j-1,:)=Savematrix1(j,:);
        elseif j==8
            obj.elK(j-1,:)=Savematrix1(j,:);
        elseif j==11
            obj.elK(j-1,:)=Savematrix1(j,:);
        end   
    end 
    clearvars Savematrix;
    clearvars Savematrix1;
%%
%         obj.elK=[12*EIz/L^3,0,0,0,0,6*EIz/L^2,12*EIz/L^3,0,0,0,0,-6*EIz/L^2;...
%         0,EA/L,0,0,0,0,0,-EA/L,0,0,0,0;...
%         0,0,-12*EIx/L^3,-6*EIx/L^2,0,0,0,0,12*EIx/L^3,-6*EIx/L^2,0,0;...
%         0,0,-6*EIx/L^2,4*EIx/L,0,0,0,0,-6*EIx/L^2,2*EIx/L,0,0;...
%         0,0,0,0,GJ/L,0,0,0,0,0,-GJ/L,0;...
%         6*EIz/L^2,0,0,0,0,4*EIz/L,6*EIz/L^2,0,0,0,0,2*EIz/L;...
%         12*EIz/L^3,0,0,0,0,6*EIz/L^2,12*EIz/L^3,0,0,0,0,6*EIz/L^2;...
%         0,-EA/L,0,0,0,0,0,EA/L,0,0,0,0;...
%         0,0,12*EIx/L^3,-6*EIx/L^2,0,0,0,0,12*EIx/L^3,-6*EIx/L^2,0,0;...
%         0,0,-6*EIx/L^2,2*EIx/L,0,0,0,0,-6*EIx/L^2,4*EIx/L,0,0;...
%         0,0,0,0,-GJ/L,0,0,0,0,0,GJ/L,0;...
%         -6*EIz/L^2,0,0,0,0,2*EIz/L,6*EIz/L^2,0,0,0,0,4*EIz/L;...
%         ];
    
    
%         Q1=[12*EIz/L^3,0,0,0,0,-6*EIz/L^2;
%         0,      140,    0,      0,      0,  0;
%         0,      0,      156,    22*L,   0,  0;
%         0,      0,      22*L,   4*L^2,  0,  0;
%         0,      0,      0,      0,      (rho*J*L)*(2/6)*420/m,  0;
%         22*L,   0,      0,      0,      0,  4*L^2;];
%     
%     Q2=[54,    0,      0,      0,      0,  -13*L;
%         0,     70,     0,      0,      0,  0;
%         0,     0,      54,     -13*L,   0,  0;
%         0,     0,      13*L,  -3*L^2, 0,  0;
%         0,     0,      0,      0,      (rho*J*L)*(1/6)*420/m,  0;
%         13*L,  0,      0,      0,      0,  -3*L^2;];
%     
%     Q3=[54,    0,      0,      0,      0,  13*L;
%         0,      70,     0,      0,      0,  0;
%         0,      0,      54,     13*L,   0,  0;
%         0,      0,      -13*L,  -3*L^2, 0,  0;
%         0,      0,      0,      0,      (rho*J*L)*(1/6)*420/m,  0;
%         -13*L,  0,      0,      0,      0,  -3*L^2;];
%     
%     Q4=[156,    0,      0,      0,      0,  -22*L;
%         0,      140,    0,      0,      0,  0;
%         0,      0,      156,    -22*L,   0,  0;
%         0,      0,      -22*L,   4*L^2,  0,  0;
%         0,      0,      0,      0,      (rho*J*L)*(2/6)*420/m,  0;
%         -22*L,   0,      0,      0,      0,  4*L^2;];
    
    
    %Calculate Element Stiffness Matrices in Global Coordinates
    obj.elKglobal=obj.T'*obj.elK*obj.T;
    
    %Ke=[12*EIz/L^3,0,0,0,0,-6*EIz/L^2,-12*EIz/L^3,0,0,0,0,-6*EIz/L^2;...
    %    0,EA/L,0,0,0,0,0,-EA/L,0,0,0,0;...
    %    0,0,12*EIx/L^3,-6*EIx/L^2,0,0,0,0,-12*EIx/L^3,-6*EIx/L^2,0,0;...
    %    0,0,-6*EIx/L^2,4*EIx/L,0,0,0,0,6*EIx/L^2,2*EIx/L,0,0;...
    %    0,0,0,0,GJ/L,0,0,0,0,0,-GJ/L,0;...
    %    -6*EIz/L^2,0,0,0,0,4*EIz/L,6*EIz/L^2,0,0,0,0,2*EIz/L;...
    %    -12*EIz/L^3,0,0,0,0,6*EIz/L^2,12*EIz/L^3,0,0,0,0,6*EIz/L^2;...
    %    0,-EA/L,0,0,0,0,0,EA/L,0,0,0,0;...
    %    0,0,-12*EIx/L^3,6*EIx/L^2,0,0,0,0,12*EIx/L^3,6*EIx/L^2,0,0;...
    %    0,0,-6*EIx/L^2,2*EIx/L,0,0,0,0,6*EIx/L^2,4*EIx/L,0,0;...
    %    0,0,0,0,-GJ/L,0,0,0,0,0,GJ/L,0;...
    %    -6*EIz/L^2,0,0,0,0,2*EIz/L,6*EIz/L^2,0,0,0,0,4*EIz/L;...
    %    ];
end
