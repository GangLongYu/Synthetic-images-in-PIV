%GENERATE generate particle images

close all; clear; clc; drawnow
fontsize = 24;
fontname = 'Times New Roman';

%% parameter setting
img_size=512;                        % image size, e.g. 512x512
rho=0.6;
rho_wall = 0.03;
particle_num=round(img_size^2*rho);  % num of particles
wall_num=round(img_size^2*rho_wall); % num of particles on the wall
z0=0.333;                            % light sheet thickness
z_move=10;                           % out-of-plane-movement when generating intensity distribution in z direction
dp=2;                                % particle diameter: unit: pixel
ddp=dp/2;                            % particle diameter variation: unit: pixel
sizey=img_size; 
sizex=img_size;
delta_x = 0.3;                       % length of chosen part of aerofoil
thick_y = sizey/sizex*delta_x;       % velocity field thickness
level = 2^16-1;                      % uint16
noise_mag = level/1000;              % noise magnitude
angle = 10;                          % rotation angle
flag = 0;                            % flag = 1 means y >= 0 and flag = 0 means y < 0
offset_ratio = 2;                    % 1m length in real world means offset_ratio pixel in image

%% load and show aerofoil and velocity field
% show aerofoil and velocity field
afoil = importdata('naca4412 aerofoil.txt','\t',1);
afoil_x = afoil.data(:,3); afoil_y = afoil.data(:,4);
vfield = importdata('NACA4412 flow field',',',1);
field_x = vfield.data(:,2); field_y = vfield.data(:,3);
field_u = vfield.data(:,4); field_v = vfield.data(:,5);
f = figure('Name','Aerofoil and Velocity'); figure(f)
plot(afoil_x,afoil_y); axis image
hold on; quiver(field_x,field_y,field_u,field_v,'g','AutoScaleFactor',0.5); hold off
xlabel('x'); ylabel('y'); title('all'); xlim([-7 15]); ylim([-7 7])
set(gcf,'position',get(0,'ScreenSize')); set(gca,'FontName',fontname,'FontSize',fontsize)
saveas(gcf,'img/all_flow_field','png');

% show part of the aerofoil and velocity field
chord = 1;
end_x = delta_x + (chord-delta_x)*rand; % randomly choose part of aerofoil
start_x = end_x - delta_x;

partfield_idx = field_x >= start_x & field_x <= end_x;
partfoil_idx = afoil_x >= start_x & afoil_x <= end_x;
if flag == 1 % the upper part of aerofoil
    partfield_idx = partfield_idx & field_y >= 0 & field_y <= thick_y;
    afoil_idx = afoil_y >= 0;
    partfoil_idx = partfoil_idx & afoil_idx;
elseif flag == 0 % the upper part of aerofoil
    thick_y = -thick_y;
    partfield_idx = partfield_idx & field_y < 0 & field_y >= thick_y;
    afoil_idx = afoil_y < 0;
    partfoil_idx = partfoil_idx & afoil_idx;
end

afoil_interp_x = afoil_x(afoil_idx);
afoil_interp_y = afoil_y(afoil_idx);
[afoil_interp_x,I] = sort(afoil_interp_x); % The grid vectors must be strictly monotonically increasing when using griddedInterpolant
afoil_interp_y = afoil_interp_y(I);
F_afoil = griddedInterpolant(afoil_interp_x,afoil_interp_y,'linear','nearest');
output_y = F_afoil([start_x, end_x]);
start_y = output_y(1); end_y = output_y(2);
if flag == 1
    afoil_xoi = [end_x,afoil_x(partfoil_idx)',start_x];
    afoil_yoi = [end_y,afoil_y(partfoil_idx)',start_y];
elseif flag == 0   
    afoil_xoi = [start_x,afoil_x(partfoil_idx)',end_x];
    afoil_yoi = [start_y,afoil_y(partfoil_idx)',end_y];
end
% part of flow field
part_x = field_x(partfield_idx); part_y = field_y(partfield_idx);
part_u = field_u(partfield_idx); part_v = field_v(partfield_idx);

min_x = start_x;  max_x = end_x;
min_y = min(0, thick_y); max_y = max(0, thick_y);
f = figure('Name','Part of Aerofoil'); figure(f)
plot(afoil_xoi,afoil_yoi,'LineWidth',3); axis image
hold on; quiver(part_x,part_y,part_u,part_v,'g','AutoScaleFactor',0.5); hold off
xlabel('x'); ylabel('y'); title('Part of Velocity Field of NACA 4412')
set(gcf,'position',get(0,'ScreenSize')); set(gca,'FontName',fontname,'FontSize',fontsize)
saveas(gcf,'img/part_flow_field','png');

%% generate particle images
disp(['Generating random artificial PIV images with ' num2str(particle_num) ' particles...'])
A=zeros(sizey,sizex); B=A;

[I0, I1, d] = I_d(particle_num,z_move,level,dp,ddp,z0);

% paritcle random position in real world
x0=min_x + rand(particle_num,1)*max_x;
y_low_limit = F_afoil(x0);
y0=y_low_limit+rand(particle_num,1).*(thick_y-y_low_limit);

% transfer real-world coordinate to image coordinate using scaling and translation
multiple = sizex / delta_x;
afoil_x_img = afoil_xoi*multiple;
x_move = min(afoil_x_img); % y need not to move
afoil_x_img = afoil_x_img-x_move;
x0 = x0*multiple-x_move;
if flag == 1
    afoil_y_img = afoil_yoi*multiple;
    y0 = y0*multiple;
elseif flag == 0
    afoil_y_img = -afoil_yoi*multiple;
    y0 = -y0*multiple;
end    
rd = -8.0 ./ d.^2;

x_interp = linspace(min_x,max_x,sizex);
y_interp = linspace(min_y,max_y,sizey);
[X_interp,Y_interp] = ndgrid(x_interp,y_interp);
[X,Y] = ndgrid(1:sizex,1:sizey);
F_u = scatteredInterpolant(part_x,part_y,part_u,'linear','nearest');
offsetx_real = F_u(X_interp,Y_interp);
F_u = griddedInterpolant(X,Y,offsetx_real,'linear','nearest'); % image coordinate
offsetx = F_u(x0,y0);
F_v = scatteredInterpolant(part_x,part_y,part_v,'linear','nearest');
offsety_real = F_v(X_interp,Y_interp);
F_v = griddedInterpolant(X,Y,offsety_real,'linear','nearest');
offsety = F_v(x0,y0);
offsetx = offsetx*offset_ratio; offsety = offsety*offset_ratio;

[xlimit1, xlimit2, ylimit1, ylimit2] = cal_extent(particle_num,x0,y0,d,sizex,sizey);                 % original image
[xlimit3, xlimit4, ylimit3, ylimit4] = cal_extent(particle_num,x0,y0,d,sizex,sizey,offsetx,offsety); % shifted image

ctr=0;
for n=1:particle_num % calculate grayscale of particle images
    ctr=ctr+1;
    if ctr==10000
        ctr=0;
        fprintf('.')
    end
    r = rd(n);
    for j=xlimit1(n):xlimit2(n) % place particles with gaussian intensity profile
        for i=ylimit1(n):ylimit2(n)
            A(i,j)=A(i,j)+I0(n)*exp(((j-x0(n))^2+(i-y0(n))^2)*r);
        end
    end
    for j=xlimit3(n):xlimit4(n)
        for i=ylimit3(n):ylimit4(n)            
            B(i,j)=B(i,j)+I1(n)*exp(((j-x0(n)+offsetx(n))^2+(i-y0(n)+offsety(n))^2)*r); 
        end
    end
end
%% add bright spots on the wall
x_wall = min_x + rand(wall_num,1)*max_x;
y_wall = F_afoil(x_wall);

% coordinate transformation
x_wall = x_wall*multiple-x_move; y_wall = y_wall*multiple;
dy = 2; % unit: pixel
y_wall = y_wall + 2*dy*rand(wall_num,1)-dy;
if flag == 0
    y_wall = -y_wall;
end 

[I0, I1, d] = I_d(wall_num,z_move,level,dp,ddp,z0);
[xlimit5, xlimit6, ylimit5, ylimit6] = cal_extent(wall_num,x_wall,y_wall,d,sizex,sizey); % original image
for n = 1:wall_num
    for j = xlimit5(n):xlimit6(n)
        for i = ylimit5(n):ylimit6(n)
            A(i,j)=A(i,j)+I0(n)*exp(((j-x_wall(n))^2+(i-y_wall(n))^2)*r);
            B(i,j)=B(i,j)+I1(n)*exp(((j-x_wall(n))^2+(i-y_wall(n))^2)*r);
        end
    end
end

% add Gaussian noise
A = A + noise_mag*randn(size(A));
B = B + noise_mag*randn(size(B));
A(A>level)=level;
B(B>level)=level;
img_fixed=uint16(A);
img_move=uint16(B);

%% test synthetic images using PIV algorithm
f = figure('Name','image1'); figure(f)
if flag == 1 % different coordinates between image and real word, flipping
    subplot(121); imshow(flip(img_fixed)); title('original image')
    subplot(122); imshow(flip(img_move)); title('shifted image')
elseif flag == 0
    subplot(121); imshow(img_fixed); title('original image')
    subplot(122); imshow(img_move); title('shifted image')
end
set(gca,'FontName',fontname,'FontSize',fontsize); set(gcf, 'position', get(0,'ScreenSize')); 
saveas(gcf,'img/synthetic_image','png');

[x,y,u,v] = PIV_test(img_move,img_fixed);
F_u = scatteredInterpolant(x0,y0,offsetx,'linear','nearest');
u_real = F_u(x,y);
F_v = scatteredInterpolant(x0,y0,offsety,'linear','nearest');
v_real = F_v(x,y);

f = figure('Name','result');figure(f);
subplot(121); imshow(img_fixed); axis image; title('original image')
hold on; quiver(x,y,u,v,'g','AutoScaleFactor',1); hold off; title('PIV result')
subplot(122); imshow(img_fixed); axis image; title('shifted image')
hold on; quiver(x,y,u_real,v_real,'g','AutoScaleFactor',1); hold off; title('real result')
set(gca,'FontName',fontname,'FontSize',fontsize); set(gcf, 'position', get(0,'ScreenSize'));
saveas(gcf,'img/result','png');

f = figure('Name','Comparison between u');figure(f);
u = u(:); u_real = u_real(:);
idx = ~isnan(u) & u_real ~= 0;
u = u(idx); u_real = u_real(idx);
plot(1:length(u),u,1:length(u_real),u_real,':');
legend('PIV','Real'); xlabel('N'); ylabel('U'); title('Comparison between u')
set(gca,'FontName',fontname,'FontSize',fontsize); set(gcf, 'position', get(0,'ScreenSize'));
saveas(gcf,'img/result_comparison','png');

% rotate effect
f = figure('Name','original image rotate'); figure(f)
img_fixed_rot = imrotate(img_fixed,angle,'loose');
img_move_rot = imrotate(img_move,angle,'loose');
subplot(131); imshow(img_fixed_rot)
subplot(132); imshow(img_move_rot)
subplot(133); imshow(img_fixed_rot)
[x,y,u,v] = PIV_test(img_move_rot,img_fixed_rot);
hold on; quiver(x,y,u,v,'g','AutoScaleFactor',1); hold off; title('PIV result')
set(gcf, 'position', get(0,'ScreenSize'))
saveas(gcf,'img/synthetic_image_rot','png');

fprintf('\n\n');



function [I0, I1, d] = I_d(num,z_move,level,dp,ddp,z0)
%I_DISTRI generate intensity distribution in z direction and diameter distribution

    % the position of Z direction
    Z0_pre=randn(num,1);  % normal distributed sheet intensity
    Z1_pre=randn(num,1);  % normal distributed sheet intensity
    Z0=Z0_pre*(z_move/200+0.5)+Z1_pre*(1-(z_move/200+0.5));
    Z1=Z1_pre*(z_move/200+0.5)+Z0_pre*(1-(z_move/200+0.5));
    I0=level*exp(-8*Z0.^2./z0^2); % the factor I0 is a function of the particle’s position Z, within the light sheet
    I0(I0>level)=level;           % uint16, 2^16-1
    I0(I0<0)=0;
    I1=level*exp(-8*Z1.^2./z0^2);
    I1(I1>level)=level;
    I1(I1<0)=0;
    
    % particle diameter distribution
    d=randn(num,1)/2; 
    d=dp+d*ddp;
    d(d<0)=0;
end



function [xlimit1, xlimit2, ylimit1, ylimit2] = cal_extent(num,x0,y0,d,sizex,sizey,offsetx,offsety)
%CAL_EXTENT calculate particle extents for images
    
    if nargin == 6
        offsetx = zeros(num,1);
        offsety = offsetx;
    end

    xlimit1=zeros(num,1);
    xlimit2=xlimit1;
    ylimit1=xlimit1;
    ylimit2=xlimit1;
    for n=1:num
        xlimit1(n)=floor(x0(n)-d(n)/2-offsetx(n)); % x min particle extent image2
        xlimit2(n)=ceil(x0(n)+d(n)/2-offsetx(n));  % x max particle extent image2
        ylimit1(n)=floor(y0(n)-d(n)/2-offsety(n)); % y min particle extent image2
        ylimit2(n)=ceil(y0(n)+d(n)/2-offsety(n));  % y max particle extent image2
    end
    xlimit1(xlimit1<1)=1;
    xlimit2(xlimit2>sizex)=sizex;
    ylimit1(ylimit1<1)=1;
    ylimit2(ylimit2>sizey)=sizey;
end



function [x,y,u,v] = PIV_test(img_fixed,img_move)
%PIV_TEST to test whether the synthetic images are correct or not

    disp('Performing PIV analysis with deforming windows and 4 passes...')
    % Standard PIV Settings
    s = cell(10,2); % To make it more readable, let's create a "settings table"
    % Parameter                         % Setting          % Options
    s{1,1}= 'Int. area 1';              s{1,2}=64;         % window size of first pass
    s{2,1}= 'Step size 1';              s{2,2}=s{1,2}/2;   % step of first pass
    s{3,1}= 'Subpix. finder';           s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                     s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                      s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';            s{6,2}=4;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';              s{7,2}=32;         % second pass window size
    s{8,1}= 'Int. area 3';              s{8,2}=16;         % third pass window size
    s{9,1}= 'Int. area 4';              s{9,2}=16;         % fourth pass window size
    s{10,1}='Window deformation';       s{10,2}='*spline'; % '*spline' is more accurate, but slower
    s{11,1}='Repeated Correlation';     s{11,2}=0;         % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
    s{12,1}='Disable Autocorrelation';  s{12,2}=0;         % 0 or 1 : Disable Autocorrelation in the first pass. 
    s{13,1}='Correlation style';        s{13,2}=0;         % 0 or 1 : Use circular correlation (0) or linear correlation (1). 
    % Standard image preprocessing settings
    p = cell(8,1);
    % Parameter                      % Setting             % Options
    p{1,1}= 'ROI';                   p{1,2}=s{5,2};        % same as in PIV settings
    p{2,1}= 'CLAHE';                 p{2,2}=1;             % 1 = enable CLAHE (contrast enhancement), 0 = disable
    p{3,1}= 'CLAHE tile number';     p{3,2}=[8 8];         % CLAHE window size
    p{4,1}= 'Highpass';              p{4,2}=0;             % 1 = enable highpass, 0 = disable
    p{5,1}= 'Highpass size';         p{5,2}=15;            % highpass size
    p{6,1}= 'Clipping';              p{6,2}=1;             % 1 = enable clipping, 0 = disable
    p{7,1}= 'Wiener';                p{7,2}=0;             % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
    p{8,1}= 'Wiener size';           p{8,2}=3;             % Wiener2 window size
    p{9,1}= 'Minimum intensity';     p{9,2}=0.0;           % Minimum intensity of input image (0 = no change) 
    p{10,1}='Maximum intensity';     p{10,2}=1.0;          % Maximum intensity on input image (1 = no change)

    disp('Performing PIV analysis with deforming windows and 4 passes...')
    img_move = preproc_PIV (img_move,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2});
    img_fixed = preproc_PIV (img_fixed,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2});
    [x,y,u,v,~] = FFT_multi (img_fixed,img_move,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2},s{11,2},s{12,2},s{13,2});
    % Remove less reliable values at the borders of the analysis
    u(:,1)=[];u(:,end)=[];u(1,:)=[];u(end,:)=[];
    v(:,1)=[];v(:,end)=[];v(1,:)=[];v(end,:)=[];
    x(:,1)=[];x(:,end)=[];x(1,:)=[];x(end,:)=[];
    y(:,1)=[];y(:,end)=[];y(1,:)=[];y(end,:)=[];   
end
