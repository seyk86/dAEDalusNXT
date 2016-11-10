%coordinates = [x,y,z] 
%size        = [sx,sy,sz]
%orientation = [alpha,beta,gamma]   0 <= alpha,beta,gamma <= 2*pi
%color       = [corner1,corner2,corner3,corner4,corner5,corner6,corner7,corner8]  0 <= b,g,r <= 1
%transparency= t                    0 <= t <= 1
%LineOnOff   = 1 (On)  or 0 (Off)

%Example cubehandles = plotcube([0,0,0],[1,1,1],[pi/6,pi/3,0],[1 1 1 1 1 1 1 1],0.5,1);

function handle = plotcube(coordinates,size,orientation,color8,transparency,LineOnOff)

    %
    %    c7---------c8
    %   /|          /|
    %  c3---------c4 |
    %  | |         | |
    %  | |         | |
    %  | c5--------c6
    %  | /         |/
    %  c1---------c2
    %


    %compute the 8 corners
    c1 = get3Dpoints(coordinates,[-size(1)/2;-size(2)/2;-size(3)/2],orientation);
    c2 = get3Dpoints(coordinates,[+size(1)/2;-size(2)/2;-size(3)/2],orientation);
    c3 = get3Dpoints(coordinates,[-size(1)/2;-size(2)/2;+size(3)/2],orientation);
    c4 = get3Dpoints(coordinates,[+size(1)/2;-size(2)/2;+size(3)/2],orientation);

    c5 = get3Dpoints(coordinates,[-size(1)/2;+size(2)/2;-size(3)/2],orientation);
    c6 = get3Dpoints(coordinates,[+size(1)/2;+size(2)/2;-size(3)/2],orientation);
    c7 = get3Dpoints(coordinates,[-size(1)/2;+size(2)/2;+size(3)/2],orientation);
    c8 = get3Dpoints(coordinates,[+size(1)/2;+size(2)/2;+size(3)/2],orientation);

    %plot the 6 sides
    hold on;
    handle(1) = fill3([c1(1) c2(1) c4(1) c3(1) c3(1) c3(1) c3(1) c3(1)],[c1(2) c2(2) c4(2) c3(2) c3(2) c3(2) c3(2) c3(2)],[c1(3) c2(3) c4(3) c3(3) c3(3) c3(3) c3(3) c3(3)],[color8(1) color8(2) color8(4) color8(3) color8(3) 0 1 color8(3)]);
    handle(2) = fill3([c2(1) c4(1) c8(1) c6(1) c6(1) c6(1) c6(1) c6(1)],[c2(2) c4(2) c8(2) c6(2) c6(2) c6(2) c6(2) c6(2)],[c2(3) c4(3) c8(3) c6(3) c6(3) c6(3) c6(3) c6(3)],[color8(2) color8(4) color8(8) color8(6) color8(6) 0 1 color8(6)]);
    handle(3) = fill3([c5(1) c6(1) c8(1) c7(1) c7(1) c7(1) c7(1) c7(1)],[c5(2) c6(2) c8(2) c7(2) c7(2) c7(2) c7(2) c7(2)],[c5(3) c6(3) c8(3) c7(3) c7(3) c7(3) c7(3) c7(3)],[color8(5) color8(6) color8(8) color8(7) color8(7) 0 1 color8(7)]);
    handle(4) = fill3([c1(1) c3(1) c7(1) c5(1) c5(1) c5(1) c5(1) c5(1)],[c1(2) c3(2) c7(2) c5(2) c5(2) c5(2) c5(2) c5(2)],[c1(3) c3(3) c7(3) c5(3) c5(3) c5(3) c5(3) c5(3)],[color8(1) color8(3) color8(7) color8(5) color8(5) 0 1 color8(5)]);
    handle(5) = fill3([c1(1) c2(1) c6(1) c5(1) c5(1) c5(1) c5(1) c5(1)],[c1(2) c2(2) c6(2) c5(2) c5(2) c5(2) c5(2) c5(2)],[c1(3) c2(3) c6(3) c5(3) c5(3) c5(3) c5(3) c5(3)],[color8(1) color8(2) color8(6) color8(5) color8(5) 0 1 color8(5)]);
    handle(6) = fill3([c3(1) c4(1) c8(1) c7(1) c7(1) c7(1) c7(1) c7(1)],[c3(2) c4(2) c8(2) c7(2) c7(2) c7(2) c7(2) c7(2)],[c3(3) c4(3) c8(3) c7(3) c7(3) c7(3) c7(3) c7(3)],[color8(3) color8(4) color8(8) color8(7) color8(7) 0 1 color8(7)]);

    %set transparency
    alpha(handle,transparency)
    %set LineStyle
    if LineOnOff == 0
        for cnt1 = 1:1:6
            set(handle(cnt1),'LineStyle','none');
        end
    end


function erg_coordinates = get3Dpoints(vec,coordinates,orientation)

    rot_mat = [...
                cos(orientation(1))*cos(orientation(2)),                  -cos(orientation(1))*sin(orientation(2)),                    sin(orientation(1));...
                cos(orientation(3))*sin(orientation(2))+sin(orientation(3))*sin(orientation(1))*cos(orientation(2)), cos(orientation(3))*cos(orientation(2))-sin(orientation(3))*sin(orientation(1))*sin(orientation(2)), -sin(orientation(3))*cos(orientation(1));...
                sin(orientation(3))*sin(orientation(2))-cos(orientation(3))*sin(orientation(1))*cos(orientation(2)), sin(orientation(3))*cos(orientation(2))+cos(orientation(3))*sin(orientation(1))*sin(orientation(2)),  cos(orientation(3))*cos(orientation(1))...
            ];
        
        a=orientation(1);
        b=orientation(2);
        c=orientation(3);
        
        Lx=[1       0       0
    0   cos(a)  -sin(a)
    0   sin(a) cos(a)];

Ly=[cos(b) 0 sin(b)
    0      1    0
    -sin(b)  0   cos(b)];

Lz=[cos(c) -sin(c)   0
    sin(c) cos(c)  0
    0           0   1];

rot_mat=Lx*Ly*Lz;

    erg_coordinates = (vec'+rot_mat*coordinates)';
