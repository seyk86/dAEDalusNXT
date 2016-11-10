%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [nasData]=read_pch(input)
%% read pch
fileID=fopen(input,'r');
nasData=[];
tline = 'dummy';
i=0;
nasLine=[];
while ischar(tline)
    tline = fgets(fileID);
    if tline(1)~='$' && ~any(tline==-1)
        if tline(1:6)~='-CONT-'
            if i~=0
                 nasData(i,:)=nasLine;
            end
            i=i+1;
            nasLine=str2double(strsplit(tline));
        else 
            nasLine=[nasLine str2double(strsplit(tline))];
        end
    end
end
nasData(i,:)=nasLine;
fclose(fileID);
nasData=nasData(:,[4:6 10:12]);

end
