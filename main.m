% Qingjun Chang @USI
%% MLC & MLC FOR DEFORMATION / MATLAB VERSION


clear,clc
close all
clear global
% clf
%


addpath(genpath([pwd '/gptoolbox']))


global v c_handle plog h_tmesh UV selectedlist h_sel coors keyflag...
    method_type icoors mlcoors me2coors mvcoors me1coors mlcoors_noS W H hbcoors
if exist([pwd '/data'], 'dir') ~= 7
    mkdir([pwd '/data']);
end
% load image
imgfile = 'images/default.jpg';
[im,~,alpha] = imread(imgfile);
[H,W]=size(im);
h_im= imshow(im);
if ~isempty(alpha)
    set(h_im,'alphadata',alpha);
end
hold on

% to select control points
if 0
    [v,v_handle] = cqj_drawPolygon();
    pause(0.1);
    delete(v_handle);
else
    load('G.mat');      % v:G v1:L v3: convex  cactus
end

plog = plot(polyshape(v'));  % control polygon
set(plog,'FaceAlpha',0.1,'LineStyle','-');
n = length(v);

selectedlist = false(1,n);

origin_v = v;

v0 = v;

v0(1,:) = v(1,:)/size(im,2);
v0(2,:) = v(2,:)/size(im,1);
[tri,xx,yy,boundary_markers] = cqj_triangle(v0',0.0001);
innerpoints = [xx,yy]';


FVC = zeros(length(xx),3);
for i = 1:length(xx)
    % Due to the rough processing here, the visualization image in the 
    % MATLAB interface is blurry, but a clear image can be obtained in 
    % MeshLab (or other visualization tools) through .obj.
    FVC(i,:) = im(round(yy(i)*size(im,1)),round(xx(i)*size(im,2)),:);
end
options.face_vertex_color = FVC;
view([0 -90])
h_tmesh = plot_mesh([xx*size(im,2),yy*size(im,1),zeros(length(xx),1)], ...
    tri(:,[2:4]),options);

set(h_tmesh,"FaceAlpha",1,'LineStyle','none');
% generate a mesh
h_sel = plot(v(1,:),v(2,:),'o','MarkerFaceColor','red');
set(h_sel,'xdata',[],'ydata',[]);
c_handle = plot(v(1,:),v(2,:),'o','MarkerFaceColor','none',...
    'MarkerSize',10);    % key release function

set(gcf,'WindowKeyPressFcn',{@keypressfcn});   % key press function
set(gcf,'WindowKeyReleaseFcn',{@keyreleasefcn});  % key release function
keyflag = false;
%
delete(h_im);
im(repmat(alpha<10,1,1,3)) = 255;
imwrite(im,'data/temp.png');

% set a mtl file
fid = fopen('data/texture.obj.mtl','w');
fprintf(fid,'# Generated by MATLAB\n');
fprintf(fid,'# Wavefront material file\n');
fprintf(fid,'newmtl material_0\n');
fprintf(fid,'Ka 0.200000 0.200000 0.200000\n');
fprintf(fid,'Kd 0.752941 0.752941 0.752941\n');
fprintf(fid,'Ks 1.000000 1.000000 1.000000\n');
fprintf(fid,'Tr 1.000000\n');
fprintf(fid,'illum 2\n');
fprintf(fid,'Ns 0.000000\n');
fprintf(fid,'map_Kd %s\n',['./temp.png']);
fclose(fid);

UV = [xx,1-yy];
writeOBJ('data/temp.obj',[xx*size(im,2),yy*size(im,1),...
    zeros(length(xx),1)],...
    tri(:,[2:4]),UV);
fid = fopen('data/temp.obj','a');
fprintf(fid,'mtllib ./texture.obj.mtl\n');
fclose(fid);

mlcoors = gbcoordinates(1,v0,innerpoints,boundary_markers);
me2coors = gbcoordinates(2,v0,innerpoints,boundary_markers);
mvcoors = gbcoordinates(3,v0,innerpoints,boundary_markers);
icoors = gbcoordinates(4,v0,innerpoints,boundary_markers,1);
me1coors = gbcoordinates(5,v0,innerpoints,boundary_markers);
mlcoors_noS = gbcoordinates(6,v0,innerpoints,boundary_markers);
tic
[b,bc] = boundary_conditions([xx yy],tri(:,[2 3 4]),v0',1:n,[],[1:n;2:n 1]');
hbcoors = kharmonic([xx yy],tri(:,[2 3 4]),b,bc)';
toc
coors = mlcoors;
method_type = 'mlc';

uicontrol('style','push',...
    'units','pix',...
    'position',[10 50 120 40],...
    'fontsize',14,...
    'string','load target',...
    'callback',{@load_cage});

uicontrol('style','push',...
    'units','pix',...
    'position',[10 10 120 40],...
    'fontsize',14,...
    'string','harmonic',...
    'callback',{@setcoors,'hbc'});

uicontrol('style','push',...
    'units','pix',...
    'position',[130 10 120 40],...
    'fontsize',14,...
    'string','mlc',...
    'callback',{@setcoors,'mlc'});

uicontrol('style','push',...
    'units','pix',...
    'position',[130 50 120 40],...
    'fontsize',14,...
    'string','mlc without S',...
    'callback',{@setcoors,'mlcNOs'});

uicontrol('style','push',...
    'units','pix',...
    'position',[250 10 120 40],...
    'fontsize',14,...
    'string','mec2',...
    'callback',{@setcoors,'mec2'});

uicontrol('style','push',...
    'units','pix',...
    'position',[250 50 120 40],...
    'fontsize',14,...
    'string','mec1',...
    'callback',{@setcoors,'mec1'});

uicontrol('style','push',...
    'units','pix',...
    'position',[370 10 120 40],...
    'fontsize',14,...
    'string','compute mvc',...
    'callback',{@setcoors,'mvc'});

uicontrol('style','push',...
    'units','pix',...
    'position',[370 50 120 40],...
    'fontsize',14,...
    'string','compute ic',...
    'callback',{@setcoors,'ic'});

uicontrol('style','slider',...
    'position',[490 50 120 40],...
    'fontsize',14,...
    'value',1/100,...
    'string','1',...
    'callback',{@recomputeIC,v0,innerpoints,boundary_markers});

uicontrol('style','push',...
    'units','pix',...
    'position',[490 10 120 40],...
    'fontsize',14,...
    'string','reset cage',...
    'callback',{@reset,origin_v});

set(gcf, 'WindowButtonDownFcn',{@Mouse_Callback,'down'} );
function recomputeIC(varargin)
global coors h_tmesh v method_type UV icoors

k = round(varargin{1}.Value*100)

icoors = gbcoordinates(4,varargin{3},varargin{4},varargin{5},k);
coors = icoors;
method_type = 'ic';


set(h_tmesh,'Vertices',[v*coors;zeros(1,size(coors,2))]');
writeOBJ(['data/temp_de_' method_type '_' num2str(k) '.obj'],h_tmesh.Vertices,...
    h_tmesh.Faces,UV);
fid = fopen(['data/temp_de_' method_type '_' num2str(k) '.obj'],'a');
fprintf(fid,'mtllib ./texture.obj.mtl\n');
fclose(fid);
end


function setcoors(varargin)
global coors h_tmesh v method_type UV icoors mlcoors me2coors mvcoors me1coors mlcoors_noS hbcoors
method_type = varargin{3}

switch method_type
    case 'mvc'
        coors = mvcoors;
    case 'ic'
        coors = icoors;
    case 'mec1'
        coors = me1coors;
    case 'mec2'
        coors = me2coors;
    case 'mlc'
        coors = mlcoors;
    case 'mlcNOs'
        coors = mlcoors_noS;
    case 'hbc'
        coors = hbcoors;
end

set(h_tmesh,'Vertices',[v*coors;zeros(1,size(coors,2))]');
writeOBJ(['data/temp_de_' method_type '.obj'],h_tmesh.Vertices,...
    h_tmesh.Faces,UV);
fid = fopen(['data/temp_de_' method_type '.obj'],'a');
fprintf(fid,'mtllib ./texture.obj.mtl\n');
fclose(fid);
end

function coors = gbcoordinates(method,cagepolygon,innerpoints,...
    boundary_markers,K)
% method:   #1: mlc     #2: mec2   #3: mvc   #4: ic    #5: mec1  #6: mlc\s
% #7: hbc
n = size(cagepolygon,2);
lambda = zeros(n,size(innerpoints,2));
if method == 1
    tic
    distances = geodis(cagepolygon,innerpoints');
    toc
end
tic
for i = 1:size(innerpoints,2)

    if boundary_markers(i)
        a = vecnorm(cagepolygon - repmat(innerpoints(:,i),1,n));
        b = zeros(n,1);
        if any(a<1e-6)         % p is one of the vertices
            b(a<1e-6) = 1;
            lambda(:,i) = b;%b;
        else     % p is on one of the edges
            newv = (cagepolygon - repmat(innerpoints(:,i),1,n))./a;
            newv = newv + newv(:,[2:end 1]);
            a = vecnorm(newv);
            [~,c] = min(a);
            d = norm(cagepolygon(:,c)-cagepolygon(:,mod(c,n)+1));
            b(c) = norm(cagepolygon(:,mod(c,n)+1)-innerpoints(:,i))/d;
            b(mod(c,n)+1) = norm(cagepolygon(:,c)-innerpoints(:,i))/d;
            lambda(:,i) = b;%b;
        end
        continue;
    end
    %     [mvc(:,i)] = mvcoordinates(innerpoints(:,i),v0);
    switch method
        case 1      % mlc
            lambda(:,i) = mlcoordinates(innerpoints(:,i),cagepolygon,...
                distances(i,:));
        case 2      % mec-2
            lambda(:,i) = mecoordinates(innerpoints(:,i),cagepolygon,2);
        case 3      % mvc
            lambda(:,i) = mvcoordinates(innerpoints(:,i),cagepolygon);
        case 4      % ic
            lambda(:,i) = icoordinates(innerpoints(:,i),cagepolygon,K);
        case 5      % mec-1
            lambda(:,i) = mecoordinates(innerpoints(:,i),cagepolygon,1);
        case 6      % mlc\s
            lambda(:,i) = mlcoordinates(innerpoints(:,i),cagepolygon,...
                ones(1,n));
    end
end
coors = lambda;
toc
end

function reset(varargin)
global v c_handle plog coors h_tmesh selectedlist h_sel
v = varargin{3};
selectedlist(selectedlist) = false;
set(c_handle,'xdata',v(1,:),'ydata',v(2,:));
set(h_sel,'xdata',v(1,selectedlist),'ydata',v(2,selectedlist));
set(plog,'Shape',polyshape(v'));
set(h_tmesh,'Vertices',[v*coors;zeros(1,size(coors,2))]');
end

function load_cage(varargin)
global v c_handle plog coors h_tmesh UV selectedlist h_sel method_type
[file,path,index] = uigetfile( ...
    {'*.mat','MAT-files (*.mat)'; ...
    '*.txt;*.xy;*.txt','text files (*.txt,*.xy,*.txt)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Select a File');
if ~index
    return;
end
cage = load([path file]);
msg = 'Arrays have incompatible sizes for this operation.';
switch index
    case 1
        assert(size(cage.v,2)==size(v,2), msg);
        v = cage.v;
    case 2
        assert(size(cage,1)==size(v,2), msg);
        v = cage';
end
selectedlist(selectedlist) = false;
set(c_handle,'xdata',v(1,:),'ydata',v(2,:));
set(h_sel,'xdata',v(1,selectedlist),'ydata',v(2,selectedlist));
set(plog,'Shape',polyshape(v'));
set(h_tmesh,'Vertices',[v*coors;zeros(1,size(coors,2))]');

writeOBJ(['data/temp_de_' method_type '.obj'],h_tmesh.Vertices,...
    h_tmesh.Faces,UV);
fid = fopen(['data/temp_de_' method_type '.obj'],'a');
fprintf(fid,'mtllib ./texture.obj.mtl\n');
fclose(fid);
end


function Mouse_Callback(~,~,action)
global v c_handle plog coors h_tmesh UV keyflag selectedlist h_sel...
    method_type H W
persistent ind xdata ydata

pos = get(gca,'CurrentPoint');
switch action
    case 'down'
        xdata = get(c_handle,'xdata');
        ydata = get(c_handle,'ydata');
        [minvalue,ind] = min(((xdata-pos(1))/W).^2+((ydata-pos(3))/H).^2);

        if minvalue<0.001
            selectedlist(ind) = ~selectedlist(ind);

            set(h_sel,'xdata',v(1,selectedlist),'ydata',v(2,selectedlist));
            set(gcf,...
                'WindowButtonMotionFcn',  {@Mouse_Callback,'move'},...
                'WindowButtonUpFcn',      {@Mouse_Callback,'up'});
        end
    case 'move'
        % move
        xdata(selectedlist) = xdata(selectedlist) + (pos(2)-xdata(ind));
        ydata(selectedlist) = ydata(selectedlist) + (pos(3)-ydata(ind));
        %         xdata(ind) = pos(2);
        %         ydata(ind) = pos(3);
        set(c_handle,'xdata',xdata,'ydata',ydata);
        v = [xdata;ydata];
        set(h_sel,'xdata',v(1,selectedlist),'ydata',v(2,selectedlist));

        set(plog,'Shape',polyshape(v'));

        set(h_tmesh,'Vertices',[v*coors;zeros(1,size(coors,2))]');

    case 'up'
        set(gcf,...
            'WindowButtonMotionFcn',  '',...
            'WindowButtonUpFcn',      '');
        if ~keyflag
            selectedlist(selectedlist) = false;
            set(h_sel,'xdata',v(1,selectedlist),'ydata',v(2,selectedlist));
        end
        writeOBJ(['data/temp_de_' method_type '.obj'],h_tmesh.Vertices,...
            h_tmesh.Faces,UV);
        fid = fopen(['data/temp_de_' method_type '.obj'],'a');
        fprintf(fid,'mtllib ./texture.obj.mtl\n');
        fclose(fid);
end
end

% key release function
function keyreleasefcn(~,evt)
global keyflag selectedlist h_sel v
switch evt.Key
    case {'control', 'shift'}
        keyflag = false;
        selectedlist(selectedlist) = false;
        set(h_sel,'xdata',v(1,selectedlist),'ydata',v(2,selectedlist));
    otherwise
        return;
end
end
function keypressfcn(~,evt)
global keyflag
switch evt.Key
    case {'control', 'shift'}
        keyflag = true;
    otherwise
        return;
end
end


function distances = geodis(polygon,points)

polygon = polygon';

n = size(polygon,1);

distances = -.1*ones(size(points,1),n);

[ins,~] = inpolygon(points(:,1),points(:,2),polygon(:,1),polygon(:,2));

innerpoints = points(ins,:);

[tri,x,y,~] = cqj_triangle(polygon,.001);
% boundary_markers = boundary_markers==1;
v = [x,y];
v(:,3) = 0;
faces = tri(:,2:4);
L = cotmatrix(v,faces);
M = massmatrix(v,faces, 'voronoi');

B = L / M * L;

g = pinv(full(B));

d = zeros(length(x),n);
for i = 1:n
    for j = 1:length(x)
        d(j,i) = g(i,i)+g(j,j)-2*g(i,j);
    end
end
d = d.^0.5;

nf = size(faces,1);
innerpoints(:,3) = 0;
innerpoints(:,4) = 0;
%A*********B
% *       *
%  *  v  *
%   *   *
%     C
AB = v(faces(:,2),:)-v(faces(:,1),:);
BC = v(faces(:,3),:)-v(faces(:,2),:);
CA = v(faces(:,1),:)-v(faces(:,3),:);

S = cross(AB,BC);

% VA = v(faces(:,1),:)-repmat;
t_points = innerpoints;
remainer = ones(size(innerpoints,1),1)==1;
re = remainer;
bary = zeros(size(innerpoints,1),3);
dis = zeros(size(innerpoints,1),n);
for i = 1:nf
    [in,~] = inpolygon(t_points(:,1),t_points(:,2),x(faces(i,:)),...
        y(faces(i,:)));

    inp = innerpoints(in & remainer,1:3);
    coors = zeros(size(inp,1),3);
    for j = 1:size(inp,1)
        VA = v(faces(i,1),:) - inp(j,:);
        VB = v(faces(i,2),:) - inp(j,:);
        VC = v(faces(i,3),:) - inp(j,:);
        coor = [cross(VB,BC(i,:));cross(VC,CA(i,:));cross(VA,AB(i,:))];
        coors(j,:) = coor(:,3)'/S(i,3);
    end

    dis(in & remainer,:) = coors*d(faces(i,:),:);
    remainer = remainer & ~in;
end

% dis = dis./sum(dis,2);


distances(ins,:) = dis;

end


% show coordinate function
% showbase(n,me1coors,tri(:,2:4),xx,yy,v0)