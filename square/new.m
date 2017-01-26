% This script assemble BIS in diamond structure and find the place for NIPA to connect

clear
clc

BIS_DISTANCE = 5; % nm
NIPA_DISTANCE = 0.154;

% Build diamond structure for BIS to repeat
% FCC A
A = [0 0 0;1 0 0;0 1 0;1 1 0];
AB =A*BIS_DISTANCE;
bis_number = size(AB,1);


inname = 'bis.gro'; % BIS structure
inname1 = 'nipa.gro'; % NIPA structure
inname2 = 'end.gro'; % end CH3
outname = 'bis3.gro'; % network structure with replicate BIS

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

% read in end molecule
end_info = importdata(inname2,delimiterIn,0);
end_molecule_name = end_info.textdata{4};
end_atom_name = end_info.textdata{5};

% store coordinate of carbon atoms that will connect to NIPA
active_atom =[];
%construct the bis pair
bis_pair = [];
%store the BIS coordinate
origin_coor = [];


%store the BIS coordinate and its active atom 
for i = 1:bis_number
    origin_coor = [origin_coor;repmat(AB(i,:),number_bis_atom,1)+one_bis_coor];
    active_atom = [active_atom;AB(i,:)+one_bis_coor(1,:);AB(i,:)+one_bis_coor(number_bis_atom,:)];
end

%build BIS carbon atom pairs
for i = 1:size(active_atom,1)
    for j = i:size(active_atom,1)
        if norm(active_atom(i,:)-active_atom(j,:)) < 0.99*BIS_DISTANCE &&  norm(active_atom(i,:)-active_atom(j,:)) > 0.7*BIS_DISTANCE
            bis_pair = [bis_pair;i,j];
        end
    end
end
% the above code just give half of the square.
% need following code to form a BIS square network.
bis_pair = [bis_pair;bis_pair'];


% ADD NIPA
% nipa molecule index
% start 10000 avoid nipa_index overlap with bis index
nipa_index = 10000;
end_index = 50000;
% store the connect pair of atoms
connect_pair = [];
connect_num = 0;
nipa_coor = [];
end_coor = [];

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
    end
    connect_pair(connect_num,2) = bis_pair(i,2);
end


% Add end NIPA and end CH3
for i = 1:size(bis_pair,1)
    nipa_between_pair = active_atom(bis_pair(i,2),:)-active_atom(bis_pair(i,1),:);
    % one end for connect chain
    if mod(bis_pair(i,1),2) == 0 
        end_bis = bis_pair(i,1)-1; % index of BIS connected to end chain
    else
        end_bis = bis_pair(i,1)+1;
    end
    nipa1 = active_atom(end_bis,:) - nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE;
    nipa2 = nipa1 - nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE;
    new_pos = nipa_pos(nipa1,nipa2,one_nipa_coor,mod(j,2));
    nipa_coor = [nipa_coor;new_pos];
    nipa_index = nipa_index+1;

    end_coor = [end_coor;nipa2 - nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE];
    end_index = end_index + 1;

    connect_num = connect_num+1;
    connect_pair(connect_num,1) = end_bis;
    connect_pair(connect_num,2) = nipa_index;

    connect_num = connect_num+1;
    connect_pair(connect_num,1) = nipa_index;
    connect_pair(connect_num,2) = end_index;

    % another ends
    if mod(bis_pair(i,2),2) == 0 
        end_bis = bis_pair(i,2)-1; % index of BIS connected to end chain
    else
        end_bis = bis_pair(i,2)+1;
    end
    nipa1 = active_atom(end_bis,:) + nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE;
    nipa2 = nipa1 + nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE;
    new_pos = nipa_pos(nipa1,nipa2,one_nipa_coor,mod(j,2));
    nipa_coor = [nipa_coor;new_pos];
    nipa_index = nipa_index+1;

    end_coor = [end_coor;nipa2 + nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE];
    end_index = end_index + 1;

    connect_num = connect_num+1;
    connect_pair(connect_num,1) = end_bis;
    connect_pair(connect_num,2) = nipa_index;

    connect_num = connect_num+1;
    connect_pair(connect_num,1) = nipa_index;
    connect_pair(connect_num,2) = end_index;
end
nipa_number = nipa_index - 10000;
end_number = end_index- 50000;




% generate network structure.gro
formatSpec = '%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n';

fileID = fopen(outname,'w');
fprintf(fileID,'%s\n',outname);
fprintf(fileID,'%5i\n',size(origin_coor,1)+size(nipa_coor,1)+size(end_coor,1));

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
for i = 1:end_number
    atom_index = atom_index+1;
    fprintf(fileID,formatSpec,i+bis_number+nipa_number,end_molecule_name, ...
                              end_atom_name,atom_index, ...
                              end_coor(i,:));
end
fprintf(fileID,'%10.5f%10.5f%10.5f\n',BIS_DISTANCE,BIS_DISTANCE,BIS_DISTANCE);
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
    elseif connect_pair(i,3) < 50000
        mol_num_tmp = connect_pair(i,3)-10000;
        connect_pair(i,3) = mol_num_tmp + bis_number;
        connect_pair(i,4) = bis_number*number_bis_atom + 8 + (mol_num_tmp-1)*number_nipa_atom;
    else 
        mol_num_tmp = connect_pair(i,3)-50000;
        connect_pair(i,3) = mol_num_tmp + bis_number + nipa_number;
        connect_pair(i,4) = mol_num_tmp + bis_number*number_bis_atom + nipa_number*number_nipa_atom;              
    end
end

connect_name = 'connect.txt';
formatSpec1 = '%5d%5d%5d%5d\n';

fileID = fopen(connect_name,'w');
fprintf(fileID,formatSpec1,connect_pair');
fclose(fileID);


% output molecule and atom connection
connect_name = 'mol.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    fprintf(fileID,'%5d',connect_pair(i,1));
    if connect_pair(i,3) <= bis_number | connect_pair(i,3) > bis_number + nipa_number
        fprintf(fileID,'%5d\n',connect_pair(i,3));
    end
end
fclose(fileID);

connect_name = 'atom.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    fprintf(fileID,'%5d%5d',connect_pair(i,2),connect_pair(i,4));
    if connect_pair(i,4) <= bis_number*number_bis_atom | ...
       connect_pair(i,4) > bis_number*number_bis_atom + nipa_number*number_nipa_atom
        fprintf(fileID,'\n');
    end
end
fclose(fileID);

% outpurt number of BIS and Connection
connect_name = 'BIS and connection.txt';
fileID = fopen(connect_name,'w');
fprintf(fileID,'%d\n%d\n',bis_number,size(bis_pair,1));
fclose(fileID);


% print some useful information
connect_name = 'information.txt';
fileID = fopen(connect_name,'w');
fprintf(fileID,['BIS_DISTANCE: %.2f nm\n' ...
                'number of BIS: %d\nnumber of NIPA: %d\n' ...
                'Box size: %.2f nm \n'],BIS_DISTANCE,bis_number,...
                nipa_number,BIS_DISTANCE);
fclose(fileID);


% generate poly file
phns = 'mol.txt';
fpn = fopen (phns, 'rt'); %open file

% poly store the number of NIPA in every chain
poly_num = 1;
%poly = zeros(100);

while feof(fpn) ~= 1 % feof is 1 when we reach the end of the file
     file = fgetl(fpn); % get one line of file
     num_of_nipa = regexp(file,'\d*\.?\d*','match');  % extract number from the line
     poly(poly_num) = size(num_of_nipa,2) - 2;
     poly_num = poly_num + 1;
end
fclose(fpn);

connect_name = 'poly.txt';
fileID = fopen(connect_name,'w');

for i = 1:length(poly)
    fprintf(fileID,'%5d\n',poly(i));
end
fclose(fileID);
