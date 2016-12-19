clear
clc
% Build diamond structure
load('A.mat')
A =[A;A+repmat([0 0 1],size(A,1),1)];
A =[A;A+repmat([0 1 0],size(A,1),1)];
A =[A;A+repmat([1 0 0],size(A,1),1)];
A = unique(A,'rows');
B= A+0.25;
for i =size(B,1):-1:1
    if max(B(i,:))>2
        B(i,:)=[];
    end
end
AB =[A;B];

NUMBER_BIS = 8;
MIN_RADIUS = 0.54;

% positions are available
avaA = 1:size(A,1);
avaB = size(A,1)+1:size(AB,1);
for j = 1:2:NUMBER_BIS
    if length(avaA)==0
        break;
    end
    newcross = avaA(randi(length(avaA)));
    lst(j) = newcross;
    [avaA,avaB] = renew_ava_list( AB,avaA,avaB,newcross,MIN_RADIUS*2);
    if length(avaB)==0
        break;
    end
    newcross = avaB(randi(length(avaB)));
    lst(j+1) = newcross;
    [avaA,avaB] = renew_ava_list( AB,avaA,avaB,newcross,MIN_RADIUS*2);
end
if length(lst)<8
    errordlg('not enough space','error')
end

[x,y,z] = sphere;
for l = 1:length(lst)
    surf(MIN_RADIUS*x+AB(lst(l),1),MIN_RADIUS*y+AB(lst(l),2),MIN_RADIUS*z+AB(lst(l),3));
    hold on
end
