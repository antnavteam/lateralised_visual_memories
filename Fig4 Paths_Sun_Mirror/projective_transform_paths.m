
function [x, y, sector, cam] = projective_transform_paths (filepath, fixedPoints_cam1, fixedPoints_cam2, display)

cd(filepath);
 
% % calib point in the real world
% fixedPoints_cam1 = [-60.5 41.5 ; -20.5 41.5 ; -60.5 0 ; -20.5 0]; 
% fixedPoints_cam2 = [20.5 41.5 ; 60.5 41.5 ; 20.5 0 ; 60.5 0]; 


%get the calib points in the movie ref
d=dir('*cam1_CALIB.csv');
Mcalib = dlmread(d(1).name, ';' , 1, 0);
movingPoints_cam1(:,1) = Mcalib(:, 1);
movingPoints_cam1(:,2) = - Mcalib(:, 2);

d=dir('*cam2_CALIB.csv');
Mcalib = dlmread(d(1).name, ';' , 1, 0);
movingPoints_cam2(:,1) = Mcalib(:, 1);
movingPoints_cam2(:,2) = - Mcalib(:, 2);

% calculate the calib matrix
tform_cam1 = fitgeotrans(movingPoints_cam1,fixedPoints_cam1,'projective');
tform_cam2 = fitgeotrans(movingPoints_cam2,fixedPoints_cam2,'projective');


% ------------get path CAM 1 ------

%get the path in the movie ref
d=dir('*cam1_SEC1.csv');
Mpath = dlmread(d(1).name, ';' , 1, 0);
ua=Mpath(:,1);
va=-Mpath(:,2);

d=dir('*cam1_SEC2.csv');
Mpath = dlmread(d(1).name, ';' , 1, 0);
ub=Mpath(:,1);
vb=-Mpath(:,2);

%checkplot
if display==1
    subplot(4,1,1);
    plot(ua,va); hold on;
    plot(ub,vb,'r'); hold on
    scatter(movingPoints_cam1(:,1), movingPoints_cam1(:,2))
    text(movingPoints_cam1(:,1), movingPoints_cam1(:,2),{'1' '2' '3' '4'})
    axis equal
    title('check camera1 path')
end

u_cam1 = [ua;ub];
v_cam1 = [va;vb];

sector_cam1=[];
sector_cam1(1:length(ua)) = 1;
sector_cam1(length(ua)+1 : length(u_cam1)) = 2;



% ------------get path CAM 2 ------

%get the path in the movie ref
d=dir('*cam2_SEC2.csv');
Mpath = dlmread(d(1).name, ';' , 1, 0);
ua=Mpath(:,1);
va=-Mpath(:,2);

d=dir('*cam2_SEC3.csv');
Mpath = dlmread(d(1).name, ';' , 1, 0);
ub=Mpath(:,1);
vb=-Mpath(:,2);


%checkplot
if display==1
    subplot(4,1,2);
    plot(ua,va); hold on;
    plot(ub,vb,'r'); hold on
    scatter(movingPoints_cam2(:,1), movingPoints_cam2(:,2))
    text(movingPoints_cam2(:,1), movingPoints_cam2(:,2),{'1' '2' '3' '4'})
    axis equal
    title('check camera2 path')
end

u_cam2 = [ua;ub];
v_cam2 = [va;vb];

sector_cam2=[];
sector_cam2(1:length(ua)) = 2;
sector_cam2(length(ua)+1 : length(u_cam2)) = 3;


% ---------- Apply the calib to get the path in the real world ref

[x1,y1] = transformPointsForward(tform_cam1,u_cam1,v_cam1);
[x2,y2] = transformPointsForward(tform_cam2,u_cam2,v_cam2);


if display==1
    subplot(4,1,3);
    plot(x1,y1); hold on;
    plot(x2,y2,'r'); hold on;
    plot([0 0], [-10 50],'k--')
    scatter(fixedPoints_cam1(:,1), fixedPoints_cam1(:,2)); hold on;
    scatter(fixedPoints_cam2(:,1), fixedPoints_cam2(:,2));
    axis equal
    title ('check overlap')
end


% ---------- Fine stick the path one after the other-------

x_error = x1(end) - x2(1) ;
y_error = y1(end) - y2(1) ;

x2_corrected = x2 + x_error;
y2_corrected = y2 + y_error;


if display==1
    subplot(4,1,4);
    plot(x1,y1); hold on;
    plot(x2_corrected,y2_corrected,'r'); hold on;
    plot([0 0], [-10 50],'k--')
    scatter(fixedPoints_cam1(:,1), fixedPoints_cam1(:,2)); hold on;
    scatter(fixedPoints_cam2(:,1), fixedPoints_cam2(:,2));
    axis equal
    title ('check overlap')
end



x=[x1 ; x2_corrected(2:end)]; % (2:end) to avoid that there is a double point at the overlap location
y=[y1 ; y2_corrected(2:end)];
sector = [sector_cam1 sector_cam2(2:end)]';
cam(1:length(x1) , 1) = 1;
cam(length(x1)+1 : length(x) ,1) = 2;




