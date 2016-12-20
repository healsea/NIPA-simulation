clear
clc
% Build diamond structure
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

B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>16
        B(i,:)=[];
    end
end
AB =[A;B];
AB =AB+repmat([-8 -8 -8],size(AB,1),1);

MIN_DIS = 10;
MAX_DIS = 20;
NUMBER_BIS = 8;
LENGTH_STRUCTURE = 1.54/(0.25*sqrt(3));

% atoms between MIN and MAX
avi = find(sqrt(sum(AB.*AB,2))>MIN_DIS/LENGTH_STRUCTURE &sqrt(sum(AB.*AB,2))<MAX_DIS/LENGTH_STRUCTURE);

atom = [0 0 0];
assign_atom = 1;
while assign_atom < NUMBER_BIS
    newatom = AB(avi(randi(length(avi))),:)+atom(assign_atom,:);
    dist = atom - repmat(newatom,size(atom,1),1);
    if length(find(sqrt(sum(dist.*dist,2))<MIN_DIS/LENGTH_STRUCTURE)) ~=0
      continue;
    end
    assign_atom = assign_atom+1;
    atom(assign_atom,:) = newatom;
end
[x,y,z] = sphere;
for l = 1:length(atom)
    surf(x/2*(MIN_DIS/LENGTH_STRUCTURE)+atom(l,1),y/2*(MIN_DIS/LENGTH_STRUCTURE)+atom(l,2),z/2*(MIN_DIS/LENGTH_STRUCTURE)+atom(l,3));
    hold on
end
