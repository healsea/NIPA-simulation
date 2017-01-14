function [ pos ] = nipa_pos( nipa1,nipa2,n_pos,index )
%   nipa_pos computes the 9 atoms position by given the 2 active carbon atoms' postion %'

% original and final active carbon bond
ori = n_pos(9,:) - n_pos(8,:);
fin = nipa2 - nipa1;
%find rotation angle

t = acos(dot(ori,fin)/norm(ori)/norm(fin));

% rotation matrix
rot_mat = makehgtform('axisrotate',cross(ori,fin),t);

% the 8th atom is the rotate center
n_pos = n_pos - repmat(n_pos(8,:),size(n_pos,1),1);
n_pos = [n_pos ones(size(n_pos,1),1)];


% new position through rotate and translation
% index rotate NIPA 180 alternatively to avoid overlap
if index == 1
    pos = transpose(rot_mat*transpose(n_pos));
else
% rotate 180
    %rot1 = makehgtform('axisrotate',fin,3.1416);
    %pos = transpose(rot_mat*rot1*transpose(n_pos));
    pos = transpose(rot_mat*transpose(n_pos));
end

pos(:,4)=[];
pos = pos + repmat(nipa1,size(n_pos,1),1);



end
