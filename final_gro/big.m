% This script assemble BIS in diamond structure and find the place for NIPA to connect

clear
clc

BIS_DISTANCE = 1.5; % nm
NIPA_DISTANCE = 0.154;
LENGTH_STRUCTURE = BIS_DISTANCE/(0.25*sqrt(3));
SIZE_OF_CELL = 2;

% Build diamond structure for BIS to repeat
% FCC A
A = [0 0 0;1 0 0;0 1 0;1 1 0;
     0 0 1;1 0 1;0 1 1;1 1 1;
     0.5 0.5 0;0.5 0.5 1;0.5 0 0.5;0.5 1 0.5;0 0.5 0.5;1 0.5 0.5];
A = [A;A+repmat([1 0 0],size(A,1),1)];
A = [A;A+repmat([0 1 0],size(A,1),1)];
A = [A;A+repmat([0 0 1],size(A,1),1)];
A = unique(A,'rows');
% FCC B
B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>SIZE_OF_CELL
        B(i,:)=[];
    end
end

AB =[A;B];


bis_number = size(AB,1);
new_bis_number = 64;

inname = 'bis.gro'; % BIS structure
inname1 = 'nipa.gro'; % NIPA structure
outname = 'bis3.gro'; % network structure with replicate BIS
outname1 = 'usefulbis.gro'; % network structure without replicate BIS

% read in original BIS coordinate
delimiterIn = ' ';
headerlinesIn = 2;
bis_info = importdata(inname,delimiterIn,headerlinesIn);
bis_info.data(end,:)=[];
bis_info.textdata(1:2,:)=[];
bis_info.textdata(end,:)=[];
number_bis_atom = size(bis_info.data,1);
%extract molecule name from cell
bis_molecule_name = bis_info.textdata{1,1}(2:end);
one_bis_coor = bis_info.data(:,2:4);

% read in original nipa coordinate
nipa_info = importdata(inname1,delimiterIn,headerlinesIn);
nipa_info.data(end,:)=[];
nipa_info.textdata(1:2,:)=[];
nipa_info.textdata(end,:)=[];
number_nipa_atom = size(nipa_info.data,1);
%extract molecule name from cell
nipa_molecule_name = nipa_info.textdata{1,1}(2:end);
one_nipa_coor = nipa_info.data(:,2:4);

% store coordinate of carbon atoms that will connect to NIPA
active_atom =[];
%construct the bis pair
bis_pair = [];
%store the BIS coordinate
origin_coor = [];
new_coor = [];


% this part construct the origin2new
% is_overlap matrix help to find the overlap molecule
is_overlap = [0 0 SIZE_OF_CELL;0 SIZE_OF_CELL 0;SIZE_OF_CELL 0 0;0 SIZE_OF_CELL SIZE_OF_CELL;...
              SIZE_OF_CELL 0 SIZE_OF_CELL; SIZE_OF_CELL SIZE_OF_CELL 0; SIZE_OF_CELL SIZE_OF_CELL SIZE_OF_CELL];
% molecule that will maintain and print
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

% % molecule that will maintain and print
% printmat = [1 9 11 13 15 16 17 18];
% % After deleting the repeated BIS molecule, each active atom will lead to the following index
% origin2new = [repmat([1;2],8,1);3;4;3;4;5;6;5;6;7;8;7;8;9;10;11;12;13;14;15;16];


AB = AB*LENGTH_STRUCTURE;
%store the BIS coordinate
for i = 1:bis_number
    origin_coor = [origin_coor;repmat(AB(i,:),number_bis_atom,1)+one_bis_coor];
    active_atom = [active_atom;AB(i,:)+one_bis_coor(1,:);AB(i,:)+one_bis_coor(number_bis_atom,:)];
    if any(i == printmat)
        new_coor = [new_coor;repmat(AB(i,:),number_bis_atom,1)+one_bis_coor];
    end
end

%build BIS carbon atom pairs
for i = 1:size(active_atom,1)
    for j = i:size(active_atom,1)
        if norm(active_atom(i,:)-active_atom(j,:)) < 0.99*BIS_DISTANCE &&  norm(active_atom(i,:)-active_atom(j,:)) > 0.7*BIS_DISTANCE
            bis_pair = [bis_pair;i,j];
        end
    end
end


% ADD NIPA
% nipa molecule index
% start 10000 avoid nipa_index overlap with bis index
nipa_index = 10000;
% store the connect pair of atoms
connect_pair = [];
new_connect_pair = [];
connect_num = 0;
nipa_coor = [];

for i = 1:size(bis_pair,1)
    % direction between BIS
    nipa_between_pair = active_atom(bis_pair(i,2),:)-active_atom(bis_pair(i,1),:);
    % number of NIPA per chain
    nipa_chain = ceil(norm(nipa_between_pair)/NIPA_DISTANCE/2);
    % rotation angle between NIPA chain and BIS direction
    theta = acos(norm(nipa_between_pair)/((nipa_chain*2+1)*NIPA_DISTANCE));
    % deviation of NIPA from BIS direction, distance and direction
    dev_dis = NIPA_DISTANCE*sin(theta);
    tmp = null(nipa_between_pair);
    dev_dir = transpose(tmp(:,1));

    connect_num = connect_num+1;
    connect_pair(connect_num,1) = bis_pair(i,1);
    new_connect_pair(connect_num,1) = origin2new(bis_pair(i,1));
    

    for j = 1: nipa_chain
        % generate a new NIPA
        nipa1 = active_atom(bis_pair(i,1),:) + (2*j - 1)*nipa_between_pair/(2*nipa_chain+1)+dev_dis*dev_dir;
        nipa2 = active_atom(bis_pair(i,1),:) + 2*j*nipa_between_pair/(2*nipa_chain+1);
        new_pos = nipa_pos(nipa1,nipa2,one_nipa_coor,mod(j,2));
        nipa_coor = [nipa_coor;new_pos];
        nipa_index = nipa_index+1;

        %store nipa1 into connect pair and nipa2 into following connect pair
        connect_pair(connect_num,2) = nipa_index;
        new_connect_pair(connect_num,2) = nipa_index;
        connect_num = connect_num+1;
        connect_pair(connect_num,1) = nipa_index;
        new_connect_pair(connect_num,1) = nipa_index;
    end
    connect_pair(connect_num,2) = bis_pair(i,2);
    new_connect_pair(connect_num,2) = origin2new(bis_pair(i,2));
end
nipa_number = size(nipa_coor,1)/number_nipa_atom;

% generate network structure.gro
formatSpec = '%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n';

fileID = fopen(outname,'w');
fprintf(fileID,'%s\n',outname);
fprintf(fileID,'%5i\n',size(origin_coor,1)+size(nipa_coor,1));

atom_index = 0;
for i = 1:bis_number
    for j = 1:number_bis_atom
        atom_index = atom_index + 1;
        atom_name = bis_info.textdata{j,2};
        fprintf(fileID,formatSpec,i,bis_molecule_name, ...
                                    atom_name,atom_index, ...
                                    origin_coor(atom_index,:));
    end
end
for i = 1:nipa_number
    for j = 1:number_nipa_atom
        atom_index = atom_index+1;
        atom_name = nipa_info.textdata{j,2};
        fprintf(fileID,formatSpec,i+bis_number,nipa_molecule_name, ...
                                  atom_name,atom_index, ...
                                  nipa_coor(atom_index-bis_number*number_bis_atom,:));
    end
end
fprintf(fileID,'%10.5f%10.5f%10.5f\n',LENGTH_STRUCTURE*SIZE_OF_CELL,LENGTH_STRUCTURE*SIZE_OF_CELL,LENGTH_STRUCTURE*SIZE_OF_CELL);
fclose(fileID);


% generate network structure.gro without replicated BIS
fileID = fopen(outname1,'w');
fprintf(fileID,'%s\n',outname1);
fprintf(fileID,'%5i\n',size(new_coor,1)+size(nipa_coor,1));

atom_index = 0;
for i = 1:new_bis_number
    for j = 1:number_bis_atom
        atom_index = atom_index + 1;
        atom_name = bis_info.textdata{j,2};
        fprintf(fileID,formatSpec,i,bis_molecule_name, ...
                                    atom_name,atom_index, ...
                                    new_coor(atom_index,:));
    end
end
for i = 1:nipa_number
    for j = 1:number_nipa_atom
        atom_index = atom_index+1;
        atom_name = nipa_info.textdata{j,2};
        fprintf(fileID,formatSpec,i+new_bis_number,nipa_molecule_name, ...
                                  atom_name,atom_index, ...
                                  nipa_coor(atom_index-new_bis_number*number_bis_atom,:));
    end
end
fprintf(fileID,'%10.5f%10.5f%10.5f\n',LENGTH_STRUCTURE*SIZE_OF_CELL,LENGTH_STRUCTURE*SIZE_OF_CELL,LENGTH_STRUCTURE*SIZE_OF_CELL);
fclose(fileID);

%output connect repair
connect_pair = [connect_pair(:,1),zeros(size(connect_pair,1),1),...
                connect_pair(:,2),zeros(size(connect_pair,1),1)];
for i = 1:size(connect_pair,1)
    if connect_pair(i,1) < 10000
        mol_num_tmp = ceil(connect_pair(i,1)/2);
        % when original mol num is odd, it is the first atom of molecule, else it is the 11st.
        connect_pair(i,2) = (mol_num_tmp-1)*number_bis_atom + 1 + (1-mod(connect_pair(i,1),2))*10;
        connect_pair(i,1) = mol_num_tmp;
    else
        mol_num_tmp = connect_pair(i,1)-10000;
        connect_pair(i,1) = mol_num_tmp + bis_number;
        connect_pair(i,2) = bis_number*number_bis_atom + 9 + (mol_num_tmp-1)*number_nipa_atom;
    end
    if connect_pair(i,3) < 10000
        mol_num_tmp = ceil(connect_pair(i,3)/2);
        % when original mol num is odd, it is the first atom of molecule, else it is the 11st.
        connect_pair(i,4) = (mol_num_tmp-1)*number_bis_atom + 1 + (1-mod(connect_pair(i,3),2))*10;
        connect_pair(i,3) = mol_num_tmp;
    else
        mol_num_tmp = connect_pair(i,3)-10000;
        connect_pair(i,3) = mol_num_tmp + bis_number;
        connect_pair(i,4) = bis_number*number_bis_atom + 8 + (mol_num_tmp-1)*number_nipa_atom;
    end
end

new_connect_pair = [new_connect_pair(:,1),zeros(size(new_connect_pair,1),1),...
                new_connect_pair(:,2),zeros(size(new_connect_pair,1),1)];
for i = 1:size(new_connect_pair,1)
    if new_connect_pair(i,1) < 10000
        mol_num_tmp = ceil(new_connect_pair(i,1)/2);
        % when original mol num is odd, it is the first atom of molecule, else it is the 11st.
        new_connect_pair(i,2) = (mol_num_tmp-1)*number_bis_atom + 1 + (1-mod(new_connect_pair(i,1),2))*10;
        new_connect_pair(i,1) = mol_num_tmp;
    else
        mol_num_tmp = new_connect_pair(i,1)-10000;
        new_connect_pair(i,1) = mol_num_tmp + new_bis_number;
        new_connect_pair(i,2) = new_bis_number*number_bis_atom + 9 + (mol_num_tmp-1)*number_nipa_atom;
    end
    if new_connect_pair(i,3) < 10000
        mol_num_tmp = ceil(new_connect_pair(i,3)/2);
        % when original mol num is odd, it is the first atom of molecule, else it is the 11st.
        new_connect_pair(i,4) = (mol_num_tmp-1)*number_bis_atom + 1 + (1-mod(new_connect_pair(i,3),2))*10;
        new_connect_pair(i,3) = mol_num_tmp;
    else
        mol_num_tmp = new_connect_pair(i,3)-10000;
        new_connect_pair(i,3) = mol_num_tmp + new_bis_number;
        new_connect_pair(i,4) = new_bis_number*number_bis_atom + 8 + (mol_num_tmp-1)*number_nipa_atom;
    end
end



connect_name = 'connect.txt';
new_connect_name = 'new_connect.txt';
formatSpec1 = '%5d%5d%5d%5d\n';

fileID = fopen(connect_name,'w');
fprintf(fileID,formatSpec1,connect_pair');
fclose(fileID);

fileID = fopen(new_connect_name,'w');
fprintf(fileID,formatSpec1,new_connect_pair');
fclose(fileID);


% output molecule and atom connection
connect_name = 'mol.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    fprintf(fileID,'%5d',connect_pair(i,1));
    if connect_pair(i,3) <= bis_number
        fprintf(fileID,'%5d\n',connect_pair(i,3));
    end
end

connect_name = 'atom.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    fprintf(fileID,'%5d%5d',connect_pair(i,2),connect_pair(i,4));
    if connect_pair(i,4) <= bis_number*number_bis_atom
        fprintf(fileID,'\n');
    end
end

connect_name = 'new_mol.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(new_connect_pair,1)
    fprintf(fileID,'%5d',new_connect_pair(i,1));
    if new_connect_pair(i,3) <= new_bis_number
        fprintf(fileID,'%5d\n',new_connect_pair(i,3));
    end
end

connect_name = 'new_atom.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(new_connect_pair,1)
    fprintf(fileID,'%5d%5d',new_connect_pair(i,2),new_connect_pair(i,4));
    if new_connect_pair(i,4) <= new_bis_number*number_bis_atom
        fprintf(fileID,'\n');
    end
end

% outpurt number of BIS and Connection
connect_name = 'BIS and connection.txt';
fileID = fopen(connect_name,'w');
fprintf(fileID,'%d\n%d\n',new_bis_number,size(bis_pair,1));