% This script assemble BIS in diamond structure and find the place for NIPA to connect

clear
clc

BIS_DISTANCE = 1.5; % nm
NIPA_DISTANCE = 0.154;
LENGTH_STRUCTURE = BIS_DISTANCE/(0.25*sqrt(3));


% Build diamond structure for BIS to repeat
% FCC A
load('A.mat')

% FCC B
B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>1
        B(i,:)=[];
    end
end
AB =[A;B];
% AB =AB+repmat([-0.5 -0.5 -0.5],size(AB,1),1);
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
connect_tmp = [];

formatSpec = '%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n';
fileID = fopen(outname,'w');
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
        %use for connect pair
        connect_tmp = [connect_tmp;II,order];
    end

    % output
    fprintf(fileID,formatSpec,II,molecule_name,atom,order,pos);
    end
end

%build BIS carbon atom pairs
% è¿™é‡Œçš„ï¿½?è·¯æ˜¯é€šè¿‡åˆ¤æ–­ç¢³åŽŸå­?é—´è·?æ?¥æ‰¾è¿žæŽ¥çš„pairï¼Œæˆ‘å?‘çŽ°æ¯?ä¸€ä¸ªåˆ†å­?èƒ½è¿žå››ä¸ªåˆ†å­?ï¼Œå…¶ä¸­æ¯?ä¸€ä¸ªåˆ†å­?çš„ä¸¤ä¸ªåŽŸå­?èƒ½è¿žä¸¤ä¸ªï¼ŒäºŽæ˜¯ï¿½?æ‹©å®ƒä»¬åˆ†åˆ«åŽ»è¿žè·?ç¦»è¾ƒå°?çš„ï¼Œï¿½?ä¸”å?‘çŽ°è¿™ä¸ªè·?ç¦»ä¸º1.15ï¿½ï¿½?.36
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

% store the connect pair of atoms
connect_pair = [];
connect_num = 0;

IIN = 0;
for i = 1:size(pair,1)
    % direction between BIS
    nipa_between_pair = active_atom(pair(i,2),:)-active_atom(pair(i,1),:);

    % number of NIPA per chain
    nipa_chain = ceil(norm(nipa_between_pair)/NIPA_DISTANCE/2);

    % rotation angle between NIPA chain and BIS direction
    theta = acos(norm(nipa_between_pair)/((nipa_chain*2+1)*NIPA_DISTANCE));

    % deviation of NIPA from BIS direction, distance and direction
    dev_dis = NIPA_DISTANCE*sin(theta);
    tmp = null(nipa_between_pair);
    dev_dir = transpose(tmp(:,1));

    connect_num = connect_num+1;
    connect_pair(connect_num,1:2) = connect_tmp(pair(i,1),:);
    for j = 1: nipa_chain
        nipa1 = active_atom(pair(i,1),:) + (2*j - 1)*nipa_between_pair/(2*nipa_chain+1)+dev_dis*dev_dir;
        nipa2 = active_atom(pair(i,1),:) + 2*j*nipa_between_pair/(2*nipa_chain+1);

        new_pos = nipa_pos(nipa1,nipa2,n_pos,mod(j,2));
        IIN = IIN+1;
        II = II+1;
        
        %store nipa1 into connect pair and nipa2 into following connect pair
        connect_pair(connect_num,3:4) = [II, 8 + (IIN-1)*number_atom1 + NUMBER_REPLICATE*number_atom];
        connect_num = connect_num+1;
        connect_pair(connect_num,1:2) = [II, 9 + (IIN-1)*number_atom1 + NUMBER_REPLICATE*number_atom];


        for k = 1:number_atom1
        % atom name
        atomN = B.textdata{k,2};
        % order of atom
        orderN = k + (IIN-1)*number_atom1 + NUMBER_REPLICATE*number_atom;

        % output
        fprintf(fileID,formatSpec,II,molecule_name1,atomN,orderN,new_pos(k,:));
        end
    end
    connect_pair(connect_num,3:4) = connect_tmp(pair(i,2),:);
end
fprintf(fileID,'%10.5f%10.5f%10.5f\n',LENGTH_STRUCTURE,LENGTH_STRUCTURE,LENGTH_STRUCTURE);
fclose(fileID);

%output connect pair
connect_name = 'connect.txt';
formatSpec1 = '%5d%5d%5d%5d\n';
fileID = fopen(connect_name,'w');
%fprintf(fileID,'%5s%5s%5s%5s\n','mol num1 ','atom num1 ','mol num2 ','atom num2 ');
fprintf(fileID,formatSpec1,connect_pair');
