clear
clc
close all

% constant. Length : nm
MIN_DIS = 1;
MAX_DIS = 2;
NUMBER_BIS = 8;
BOX_SIZE =4;
LENGTH_STRUCTURE = 0.154/(0.25*sqrt(3));


% Build diamond structure
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
A =[A;A+repmat([0 0 4],size(A,1),1)];
A =[A;A+repmat([0 4 0],size(A,1),1)];
A =[A;A+repmat([4 0 0],size(A,1),1)];
A = unique(A,'rows');
A =[A;A+repmat([0 0 8],size(A,1),1)];
A =[A;A+repmat([0 8 0],size(A,1),1)];
A =[A;A+repmat([8 0 0],size(A,1),1)];
A = unique(A,'rows');
% FCC B
B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>16
        B(i,:)=[];
    end
end
%AB =[A;B];
%AB =AB+repmat([-8 -8 -8],size(AB,1),1);
B = B+repmat([-8 -8 -8],size(B,1),1);
B = B*LENGTH_STRUCTURE;

% atoms between MIN and MAX in FCC_B since we choose [0 0 0] as reference point which is in FCC_A
avi = find(sqrt(sum(B.*B,2))>MIN_DIS &sqrt(sum(B.*B,2))<MAX_DIS);

% initial atom at the center of the box
atom = ones(1,3)*BOX_SIZE/2;
assign_atom = 1;
while assign_atom < NUMBER_BIS
    newatom = B(avi(randi(length(avi))),:)+atom(assign_atom,:);

        % new atom must be in the box
    if length(find(newatom<0)) ~= 0 || length(find(newatom>BOX_SIZE)) ~=0
        continue;
    end
    
    dist = atom - repmat(newatom,size(atom,1),1);
    % new atom should not be too close to previous atoms and their mirrors
    if length(find(sqrt(sum(dist.*dist,2))<MIN_DIS)) ~=0 || length(find(sqrt(sum(dist.*dist,2))>BOX_SIZE-MIN_DIS)) ~=0
        continue;
    end

    assign_atom = assign_atom+1;
    atom(assign_atom,:) = newatom;
end
scatter3(atom(:,1),atom(:,2),atom(:,3))

%[x,y,z] = sphere;
% for l = 1:length(atom)
%     surf(x/2*(MIN_DIS/LENGTH_STRUCTURE)+atom(l,1),y/2*(MIN_DIS/LENGTH_STRUCTURE)+atom(l,2),z/2*(MIN_DIS/LENGTH_STRUCTURE)+atom(l,3));
%     hold on
% end
