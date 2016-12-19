clear
clc
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
NUMBER_BIS = 8;
MIN_DISTANCE = 0.5;
crosslink =[];
for j = 1:2:NUMBER_BIS
    newcross = A(randi(size(A,1)),:);
    k = 1;
    while k <= size(crosslink,1)
       if norm(newcross-crosslink(k,:))<MIN_DISTANCE
          newcross = A(randi(size(A,1)),:);
          k = 1;
          continue;
       end
       k = k+1;
    end
    crosslink(j,:) = newcross;
    newcross = B(randi(size(B,1)),:);
    k = 1;
    while k <= size(crosslink,1)
       if norm(newcross-crosslink(k,:))<MIN_DISTANCE
          newcross = B(randi(size(B,1)),:);
          k = 1;
          continue;
       end
       k = k+1;
    end
    crosslink(j+1,:) = newcross;
end

crosslink

[x,y,z] = sphere;
for l = 1:length(crosslink)
    surf(MIN_DISTANCE*x+crosslink(l,1),MIN_DISTANCE*y+crosslink(l,2),MIN_DISTANCE*z+crosslink(l,3));
    hold on
end
