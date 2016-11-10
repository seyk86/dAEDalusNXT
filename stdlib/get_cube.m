
function [c1 c2 c3 c4 c5 c6 c7 c8] = get_cube(coordinates,size,orientation)

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

rot_mat=Lx*Lz*Ly;

    erg_coordinates = (vec'+rot_mat*coordinates)';
