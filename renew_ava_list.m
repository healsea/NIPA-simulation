function [ newavaA,newavaB ] = renew_ava_list( AB,avaA,avaB,newcross,MIN_DISTANCE)
%renew_ava_list renew the available list when a new atom is inserted
newavaA = avaA;
newavaB = avaB;
for i = 1:length(avaA)
    if (norm(AB(avaA(i),:)-AB(newcross,:))<MIN_DISTANCE)
        newavaA(i)=0;
    end
end
for i = 1:length(avaB)
    if (norm(AB(avaB(i),:)-AB(newcross,:))<MIN_DISTANCE)
        newavaB(i)=0;
    end
end
newavaA(newavaA==0)=[];
newavaB(newavaB==0)=[];
end
