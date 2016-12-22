clear
clc
close all

% constant Length : nm
NUMBER_BIS = 8;
BOX_SIZE =4;
BIS_LENGTH = 0.591;
CC_BOND =0.154;

% diamond direction
dir = [0,0,1;0,sqrt(8)/3,-1/3;sqrt(6)/3,-sqrt(2)/3,-1/3;-sqrt(6)/3,-sqrt(2)/3,-1/3];
dir = [dir;-dir]*CC_BOND;

bis_pos = [BOX_SIZE/3 BOX_SIZE/3 BOX_SIZE/3;2*BOX_SIZE/3 BOX_SIZE/3 BOX_SIZE/3;BOX_SIZE/3 2*BOX_SIZE/3 BOX_SIZE/3;BOX_SIZE/3 BOX_SIZE/3 2*BOX_SIZE/3;2*BOX_SIZE/3 2*BOX_SIZE/3 BOX_SIZE/3;2*BOX_SIZE/3 BOX_SIZE/3 2*BOX_SIZE/3;BOX_SIZE/3 2*BOX_SIZE/3 2*BOX_SIZE/3;2*BOX_SIZE/3 2*BOX_SIZE/3 2*BOX_SIZE/3];
scatter3(bis_pos(:,1),bis_pos(:,2),bis_pos(:,3),'bo');
xlim([0,BOX_SIZE]);
ylim([0,BOX_SIZE]);
zlim([0,BOX_SIZE]);
% random direction of BIS

% azimuth = 360*rand(NUMBER_BIS,1);
% polar = 180*rand(NUMBER_BIS,1)-90;
% bis_pos1 = bis_pos +[BIS_LENGTH*sin(polar).*cos(azimuth),BIS_LENGTH*sin(polar).*sin(azimuth),BIS_LENGTH*cos(polar)];
% hold on
% scatter3(bis_pos1(:,1),bis_pos1(:,2),bis_pos1(:,3),'ro');

% positions that are already occupied
occupied = [];

for i = 1:NUMBER_BIS
    % flag for jump two for loop
    flag = 0;

    azimuth = 360*rand;
    polar = 180*rand-90;
    % the other side of BIS
    bis_pos1(i,:) = bis_pos(i,:) +[BIS_LENGTH*sin(polar).*cos(azimuth),BIS_LENGTH*sin(polar).*sin(azimuth),BIS_LENGTH*cos(polar)];

    % use midpoint to check whether BIS segment is too close to other BIS or NIPA
    new_occ = [linspace(bis_pos(i,1),bis_pos1(i,1),10);linspace(bis_pos(i,2),bis_pos1(i,2),10);linspace(bis_pos(i,3),bis_pos1(i,3),10)]';%'
    for j = 1:size(occupied,1)
        dist = new_occ - repmat(occupied(j,:),size(new_occ,1),1);
        if length(find(sqrt(sum(dist.*dist,2))<CC_BOND)) ~=0 || length(find(sqrt(sum(dist.*dist,2))>BOX_SIZE-CC_BOND)) ~=0
            flag = 1;
            break;
        end
    end
    if flag == 1
        continue;
    end
    occupied = [occupied;new_occ];
end

bis_p =[bis_pos;bis_pos1];
hold on
scatter3(bis_pos1(:,1),bis_pos1(:,2),bis_pos1(:,3),'ro');

% start with an BIS
ini = randi(NUMBER_BIS*2);
start = bis_p(ini,:);
% previous to ensure the atom goes the right direction. IF here we choose
% previous as the atom connected to the knot, it would be better
previous = start;
NIPA = start;
while true
    next = start + dir(randi(8),:);
    % ensure the atom goes the right direction
    if (norm(previous - start) ~=0) & (norm(next-previous)>1.9*CC_BOND || norm(next-previous)<1.4*CC_BOND) 
        continue;
    end

    % new position may exceed box
    next = mod(next , BOX_SIZE);

    % avoid place: include previous NIPA(include start BIS) and BIS midpoint
    avoid = [NIPA;setdiff(occupied,bis_p,'rows')];
    % distance between new and avoid
    dist1 = avoid - repmat(next,size(avoid,1),1);
    if length(find(sqrt(sum(dist1.*dist1,2))<CC_BOND)) ~=0 || length(find(sqrt(sum(dist1.*dist1,2))>BOX_SIZE-CC_BOND)) ~=0
        continue;
    end

    % terminal place: ALL BIS except start BIS
    terminal = bis_p;
    terminal(ini,:)=[];
    dist2 = terminal - repmat(next,size(terminal,1),1);
    if length(find(sqrt(sum(dist2.*dist2,2)) <= CC_BOND)) ~=0 || length(find(sqrt(sum(dist2.*dist2,2)) >= BOX_SIZE-CC_BOND)) ~=0
        disp('a chain has been formed! yeah!')
        break;
    end

    % generate new NIPA
    draw = [start',next'];
    line(draw(1,:),draw(2,:),draw(3,:));
    pause(0.3);
    NIPA =[NIPA;next];
    previous = start;
    start = next;

end
