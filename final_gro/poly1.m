phns = 'new_mol.txt';
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
