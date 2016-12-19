clear
clc
load('A.mat')
A =[A;A+repmat([0 0 1],length(A),1)];
A =[A;A+repmat([0 1 0],length(A),1)];
A =[A;A+repmat([1 0 0],length(A),1)];
A = unique(A,'rows')
B= A+0.25;
for i =length(B):-1:1
    if max(B(i,:))>2
        B(i,:)=[];
    end
end
NUMBER_BIS = 8;
for j = 1:2:NUMBER_BIS
    crosslink(j,:) = A(randi(length(A)),:);
    crosslink(j+1,:) = B(randi(length(B)),:);
end
crosslink





% plot3(A(:,1),A(:,2),A(:,3),'o')
