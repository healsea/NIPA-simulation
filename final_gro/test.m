clear;
SIZE_OF_CELL = 1;

A = [0 0 0;1 0 0;0 1 0;1 1 0;
     0 0 1;1 0 1;0 1 1;1 1 1;
     0.5 0.5 0;0.5 0.5 1;0.5 0 0.5;0.5 1 0.5;0 0.5 0.5;1 0.5 0.5];

% FCC B
B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>SIZE_OF_CELL
        B(i,:)=[];
    end
end

AB =[A;B];
bis_number = 18;
% is_overlap matrix help to find the overlap molecule
is_overlap = [0 0 SIZE_OF_CELL;0 SIZE_OF_CELL 0;SIZE_OF_CELL 0 0;0 SIZE_OF_CELL SIZE_OF_CELL;...
              SIZE_OF_CELL 0 SIZE_OF_CELL; SIZE_OF_CELL SIZE_OF_CELL 0; SIZE_OF_CELL SIZE_OF_CELL SIZE_OF_CELL];
% molecule that will maintain and print
%printmat = [1 9 11 13 15 16 17 18];
printmat = [];
% After deleting the repeated BIS molecule, each active atom will lead to the following index
origin2new = zeros(bis_number*2,1);
useful_mol = 1;
% the flag judge whether the new atom is overlapped with previous atoms
flag = 0;
for i = 1:bis_number
    flag = 0;
    for j = 1:i-1
        if ismember(abs(AB(i,:)-AB(j,:)),is_overlap,'rows')
            origin2new(2*i-1:2*i) = [origin2new(2*j-1:2*j)];
            flag = 1;
            break
        end
    end
    if flag == 0
        printmat = [printmat; i];
        origin2new(2*i-1:2*i) = [2*useful_mol-1;2*useful_mol];
        useful_mol = useful_mol + 1;
    end
end
