% This script assemble BIS in diamond structure and find the place for NIPA to connect

clear
clc

BIS_DISTANCE = 1.5; % nm
NIPA_DISTANCE = 0.154;
LENGTH_STRUCTURE = BIS_DISTANCE/(0.25*sqrt(3));


% Build diamond structure for BIS to repeat
% FCC A
load('A.mat')
A =[A;A+repmat([0 0 1],size(A,1),1)];
A =[A;A+repmat([0 1 0],size(A,1),1)];
A =[A;A+repmat([1 0 0],size(A,1),1)];
A = unique(A,'rows');
A =[A;A+repmat([0 0 2],size(A,1),1)];
A =[A;A+repmat([0 2 0],size(A,1),1)];
A =[A;A+repmat([2 0 0],size(A,1),1)];
A = unique(A,'rows');
% FCC B
B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>4
        B(i,:)=[];
    end
end
AB =[A;B];
AB =AB+repmat([-2 -2 -2],size(AB,1),1);
AB = AB*LENGTH_STRUCTURE;
NUMBER_REPLICATE = size(AB,1);

inname = 'bis.gro';
inname1 = 'nipa.gro';
outname = 'bis2.gro';

% read in original bis coordinate
delimiterIn = ' ';
headerlinesIn = 2;
A = importdata(inname,delimiterIn,headerlinesIn);
A.data(end,:)=[];
A.textdata(1:2,:)=[];
A.textdata(end,:)=[];
number_atom = size(A.data,1);
%extract molecule name from cell
molecule_name = A.textdata{1,1}(2:end);

% read in original nipa coordinate
B = importdata(inname1,delimiterIn,headerlinesIn);
B.data(end,:)=[];
B.textdata(1:2,:)=[];
B.textdata(end,:)=[];
number_atom1 = size(B.data,1);
%extract molecule name from cell
molecule_name1 = B.textdata{1,1}(2:end);

% store coordinate of carbon atoms that will connect to NIPA
active_atom =[];


formatSpec = '%5i%s%5s%8i%8.3f%8.3f%8.3f\n';
fileID = fopen(outname,'a');
fprintf(fileID,'%s\n',outname);
% fprintf(fileID,'%5i\n',NUMBER_REPLICATE*number_atom);
for i = 1:NUMBER_REPLICATE
    for j = 1:number_atom
    II = i;
    % atom name
    atom = A.textdata{j,2};
    % order of atom
    order = A.data(j,1) + (i-1)*number_atom;
    %position of atom
    pos = [A.data(j,2) + AB(i,1)  ; A.data(j,3) + AB(i,2) ;A.data(j,4) + AB(i,3)];

    % store active carbon atoms
    if j == 1 || j == number_atom
        active_atom =[active_atom;pos']; %'
    end

    % output
    fprintf(fileID,formatSpec,II,molecule_name,atom,order,pos);
    end
end

fclose(fileID);


%build BIS carbon atom pairs
%% 这里的�?路是通过判断碳原子间距来找连接的pair，我发现每一个分子能连四个分子，其中每一个分子的两个原子能连两个，于是�?择它们分别去连距离较小的，�?且发现这个距离为1.15��?.36
pair = [];
for i = 1:size(active_atom,1)
    for j = i:size(active_atom,1)
        if norm(active_atom(i,:)-active_atom(j,:)) < 0.99*BIS_DISTANCE &&  norm(active_atom(i,:)-active_atom(j,:)) > 0.7*BIS_DISTANCE
        %if norm(active_atom(i,:)-active_atom(j,:)) == BIS_DISTANCE
            pair = [pair;i,j];
        end
    end
end

% put NIPA
n_pos = B.data(:,2:end);

IIN = 0;
fileID = fopen(outname,'a');
for i = 1:size(pair,1)
    nipa_between_pair = active_atom(pair(i,2),:)-active_atom(pair(i,1),:);
    for j = 1: floor(norm(nipa_between_pair)/NIPA_DISTANCE/2)
        nipa1 = active_atom(pair(i,1),:) + (2*j - 1)*nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE;
        nipa2 = nipa1 + nipa_between_pair/norm(nipa_between_pair)*NIPA_DISTANCE;
        new_pos = nipa_pos(nipa1,nipa2,n_pos);
        IIN = IIN+1;
        for k = 1:number_atom1
        % atom name
        atomN = B.textdata{k,2};
        % order of atom
        orderN = B.data(k,1) + (IIN-1)*number_atom1 + NUMBER_REPLICATE*number_atom;

        % output
        fprintf(fileID,formatSpec,IIN,molecule_name1,atomN,orderN,new_pos(k,:));
        end
    end
end
fprintf(fileID,'%10.5f%10.5f%10.5f\n',0,0,0);
fclose(fileID);
