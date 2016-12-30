% output molecule and atom connection
connect_name = 'mol1.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    fprintf(fileID,'%5d',connect_pair(i,1));
    if connect_pair(i,3) <= USEFUL_BIS
        fprintf(fileID,'%5d\n',connect_pair(i,3));
    end
end


connect_name = 'atom1.txt';
fileID = fopen(connect_name,'w');

for i = 1:size(connect_pair,1)
    fprintf(fileID,'%5d%5d',connect_pair(i,2),connect_pair(i,4));
    if connect_pair(i,4) <= USEFUL_BIS*number_atom
        fprintf(fileID,'\n');
    end
end
