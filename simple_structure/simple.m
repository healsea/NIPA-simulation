% This script assemble BIS in diamond structure and find the place for NIPA to connect

clear
clc

NIPA_DISTANCE = 0.154;
CHAIN_NUM = 10;
% the active atom of BIS on the top will contact to the upper chain 
DIRECTION = [-0.5,0.5,0.5;0.5,-0.5,0.5;-0.5,-0.5,-0.5;0.5,0.5,-0.5];
% one chain is CHAIN_NUM nipa and an end atom
LENGTH_STRUCTURE = (2*CHAIN_NUM+1)*NIPA_DISTANCE*2/sqrt(3);

inname = 'bis.gro'; % BIS structure
inname1 = 'nipa.gro'; % NIPA structure
inname2 = 'end.gro'; % end molecule
outname = 'simple.gro'; % network structure with replicate BIS

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


%store the BIS and NIPA coordinate
coor = one_bis_coor;

% ADD NIPA
% nipa molecule index
nipa_index = 10000;
end_index  = 50000;
% store the connect pair of atoms
connect_pair = [];
connect_num = 0;

for i = 1:size(DIRECTION,1)
	if i <= 2 
		bis_coor = one_bis_coor(1,:);
		connect_num = connect_num+1;
		connect_pair(connect_num,1) = 1;
	else
		bis_coor = one_bis_coor(11,:);
		connect_num = connect_num+1;
		connect_pair(connect_num,1) = 11;
	end
	for j = 1:CHAIN_NUM
		% one active part of bis is connected to 2 NIPA chain

		nipa1 = bis_coor + (2*j-1)*NIPA_DISTANCE*DIRECTION(i,:)/norm(DIRECTION(i,:));
		nipa2 = nipa1 + NIPA_DISTANCE*DIRECTION(i,:)/norm(DIRECTION(i,:));
		new_pos = nipa_pos(nipa1,nipa2,one_nipa_coor,mod(j,2));
		coor = [coor;new_pos];
		nipa_index = nipa_index + 1;

        connect_pair(connect_num,2) = nipa_index;
        connect_num = connect_num+1;
        connect_pair(connect_num,1) = nipa_index;
    end
    % add one end atom CH3
    end_coor = nipa2 + NIPA_DISTANCE*DIRECTION(i,:)/norm(DIRECTION(i,:));
    coor = [coor;end_coor];
    end_index = end_index + 1;
    connect_pair(connect_num,2) = end_index;   
end


% generate network structure.gro
formatSpec = '%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n';

fileID = fopen(outname,'w');
fprintf(fileID,'%s\n',outname);
fprintf(fileID,'%5i\n',size(coor,1));

for i = 1:number_bis_atom
    atom_name = bis_info.textdata{i,2};
    fprintf(fileID,formatSpec,1,bis_molecule_name, ...
            atom_name,i,coor(i,:));
end
for i = 1:size(DIRECTION,1)
    for j = 1:CHAIN_NUM
        for k = 1:number_nipa_atom
            atom_name = nipa_info.textdata{k,2};
            atom_index = number_bis_atom+(i-1)*CHAIN_NUM*number_nipa_atom+i-1+(j-1)*number_nipa_atom+k;
            fprintf(fileID,formatSpec,1+(i-1)*(CHAIN_NUM+1)+j,nipa_molecule_name, ...
                                      atom_name,atom_index,coor(atom_index,:));
        end
    end
    end_ind = number_bis_atom + i*CHAIN_NUM*number_nipa_atom + i;
    fprintf(fileID,formatSpec,1+i*CHAIN_NUM+i,end_molecule_name, ...
            end_atom_name,end_ind,coor(end_ind,:));
end


fprintf(fileID,'%10.5f%10.5f%10.5f\n',LENGTH_STRUCTURE,LENGTH_STRUCTURE,LENGTH_STRUCTURE);
fclose(fileID);


%output connect repair
connect_pair = [connect_pair(:,1),zeros(size(connect_pair,1),1),...
                connect_pair(:,2),zeros(size(connect_pair,1),1)];
for i = 1:size(connect_pair,1)
    if connect_pair(i,1) < 10000
        connect_pair(i,2) = connect_pair(i,1);
        connect_pair(i,1) = 1;
    else
        mol_num_tmp = connect_pair(i,1)-10000;
        connect_pair(i,1) = mol_num_tmp + 1 + ceil(mol_num_tmp/CHAIN_NUM)-1;
        connect_pair(i,2) = number_bis_atom + 9 + (mol_num_tmp-1)*number_nipa_atom + ceil(mol_num_tmp/CHAIN_NUM)-1;
    end 
    if connect_pair(i,3) < 10000
        connect_pair(i,4) = connect_pair(i,3);
        connect_pair(i,3) = 1;
    elseif (connect_pair(i,3) < 50000)        
        mol_num_tmp = connect_pair(i,3)-10000;
        connect_pair(i,3) = mol_num_tmp + 1 + ceil(mol_num_tmp/CHAIN_NUM)-1;
        connect_pair(i,4) = number_bis_atom + 8 + (mol_num_tmp-1)*number_nipa_atom + ceil(mol_num_tmp/CHAIN_NUM)-1;
    else
        mol_num_tmp = connect_pair(i,3)-50000;
        connect_pair(i,3) = mol_num_tmp*(CHAIN_NUM+1) + 1;
        connect_pair(i,4) = number_bis_atom + mol_num_tmp*CHAIN_NUM*number_nipa_atom + mol_num_tmp;        
    end
end

% generate connect file
connect_name = 'connect.txt';
formatSpec1 = '%5d%5d%5d%5d\n';

fileID = fopen(connect_name,'w');
fprintf(fileID,formatSpec1,connect_pair');
fclose(fileID);



% output molecule and atom connection
connect_name = 'mol.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    if connect_pair(i,1) <= 1 && i~=1  % avoid \n print at the first line
        fprintf(fileID,'%5d',connect_pair(i-1,3)); % print the end molecules
        fprintf(fileID,'\n');
    end
    fprintf(fileID,'%5d',connect_pair(i,1));
end
fprintf(fileID,'%5d',connect_pair(i,3)); % print the last end molecule
fclose(fileID);

connect_name = 'atom.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
	if connect_pair(i,2) <= number_bis_atom && i~=1  % avoid \n print at the first line
        fprintf(fileID,'\n');
    end
    fprintf(fileID,'%5d%5d',connect_pair(i,2),connect_pair(i,4));
end
fclose(fileID);


% print some useful information
connect_name = 'information.txt';
fileID = fopen(connect_name,'w');
fprintf(fileID,['number of BIS: %d\nnumber of NIPA: %d\n' ...
                'number of end: %d\nBox size: %.2f nm \n'],...
                1,size(DIRECTION,1)*CHAIN_NUM,size(DIRECTION,1),LENGTH_STRUCTURE);
fclose(fileID);


% Print number of NIPA in every chain
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
