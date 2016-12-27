clear
clc

NUMBER_REPLICATE = 3;
DISTANCE = 1;
inname = '001.txt';
outname = 'nipam.gro';
delimiterIn = ' ';
headerlinesIn = 2;
A = importdata(inname,delimiterIn,headerlinesIn);
A.data(end,:)=[];
A.textdata(1:2,:)=[];
A.textdata(end,:)=[];
number_atom = size(A.data,1);
%extract molecule name from cell
molecule_name = A.textdata{1,1}(2:end);

formatSpec = '%5i%5s%5s%5i%8.3f%8.3f%8.3f\n';
fileID = fopen(outname,'a');
fprintf(fileID,'nipam\n');
fprintf(fileID,'%5i\n',NUMBER_REPLICATE^2*number_atom);
for i = 1:NUMBER_REPLICATE^2
    for j = 1:number_atom
    II = i;
    % atom name
    atom = A.textdata{j,2};
    % order of atom
    order = A.data(j,1) + (i-1)*number_atom;
    %position of atom  THIS IS JUST FOR TWO DIMENSION
    pos = [A.data(j,2) + DISTANCE*(mod(i,NUMBER_REPLICATE)-1); A.data(j,3) + DISTANCE*(floor(i/NUMBER_REPLICATE));A.data(j,4)];

    % output
    fprintf(fileID,formatSpec,II,molecule_name,atom,order,pos);
    end
end
fprintf(fileID,'%10.5f%10.5f%10.5f\n',0,0,0);
fclose(fileID);
