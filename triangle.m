function [tri,x,y,boundary_markers] = triangle(polygon,s)
% CopyRight:  Qingjun Chang @USI
if nargin == 1
    s = 0.0001;
end
s = vpa(s);
s = arrayfun(@char, s, 'uniform', 0);
if ispc
    guid = uuidgen('mex');
else
    guid = '0000000000000';
end
filename = ['temp_' guid];
fid = fopen([filename '.poly'],'w');
fprintf(fid,['# Generated by MATLAB data_' filename '.poly\n']);
fprintf(fid,'%d %d %d %d\n',size(polygon,1),2,0,1);
for i = 1:size(polygon,1)
    fprintf(fid,'%d %f %f %d\n',i,polygon(i,1),polygon(i,2),1);
end
fprintf(fid,'%d %d\n',size(polygon,1),1);
ia = 1:size(polygon,1);
for i = 1:size(polygon,1)
    fprintf(fid,'%d %d %d %d\n',i,ia(mod(i-1,size(ia,2))+1),ia(mod(i,size(ia,2))+1),1);
end
fprintf(fid,'%d\n',0);
fclose(fid);
if ispc
    [status,~] = dos(['triangle_win -pqa' s{1} 'Dg ' filename '.poly']);
elseif ismac
    [status,~] = system(['./triangle_mos -pqa' s{1} 'Dg ' filename '.poly']);
elseif isunix
    [status,~] = system(['triangle_unix -pqa' s{1} 'Dg ' filename '.poly']);
else
    disp("Error");
end
if status ~=0
    delete([filename,'*']);
    return;
end

basename = [filename '.1'];
nodefile=[basename '.node'];
elefile=[basename '.ele'];
fid_nod = fopen(nodefile);
[nnode] = fscanf(fid_nod,'%i',[1 4]);
data = fscanf(fid_nod,'%f',[4 nnode(1)])';
x=data(:,2); y=data(:,3); boundary_markers = data(:,4);
fid_ele = fopen(elefile);
[nelem] = fscanf(fid_ele,'%i',[1 3]);
tri = fscanf(fid_ele,'%i',[4 nelem(1)])';
fclose(fid_nod);
fclose(fid_ele);

delete([filename,'*']);

end
