%% PROJECT 2 VELOCITY ESTIMATION
close all;
clear all;
clc;
%Change this for both dataset 1 and dataset 4. Do not use dataset 9.
datasetNum = 4;
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
%% INITIALIZE CAMERA MATRIX AND OTHER NEEDED INFORMATION
K_Mat=[311.0520, 0, 201.8724;0,311.3885, 113.6210;0, 0, 1]; 
avg_time_diff=0.02;
pointTracker = vision.PointTracker;
CTB = [-0.04;0;-0.03];
H=[0.7071, -0.7071, 0.0000, -0.0400;-0.7071, -0.7071, 0.0000, 0;0.0000, 0.0000, -1.0000, -0.0300;0, 0, 0, 1.0000];
estimated_V=zeros(6,length(sampledData)); EstimatedV=zeros(6,length(sampledData)); Estimatedv=zeros(6,length(sampledData));
for n = 2:length(sampledData)
    %% Initalize Loop load images
imgnow=sampledData(n).img;
imgprev=sampledData(n-1).img;
    %% Detect good points 
    points_1 = detectFASTFeatures(imgnow);
    points_1=selectStrongest(points_1,50);
    %% Initalize the tracker to the last frame.
    initialize(pointTracker,points_1.Location,imgprev);
    %% Find the location of the next points;
    [points_2,point_validity] = pointTracker(imgnow);
    pointTracker.release();
    %% Calculate Velocity
    prev_points=points_1.Location; new_points= points_2;
    new_points=[new_points,ones(length(new_points),1)]';
    new1_points=inv(K_Mat)*new_points;
    new2_points=[new1_points(1,:);new1_points(2,:);ones(1,length(new1_points))];
    prev_points=[prev_points,ones(length(prev_points),1)]';
    prev_points=inv(K_Mat)*prev_points;
    delta_vel=(new1_points-prev_points)./avg_time_diff;
    delta_vel=delta_vel(1:2,:);
    %% Calculate Height
    [position,orientation,R_c2w]=estimatePose(sampledData,n);
    T=CTB - R_c2w.'*position;
    R=R_c2w.';
    h=[R(:,1:2),T];
    A=h\prev_points; 
    Z=[];
    for i=1:length(prev_points)
        Z=[Z, 1/A(3,i)];
    end
    %% RANSAC
    e=0.4;                  %Change e between 0.3 to 0.8
    Vel=velocityRANSAC(delta_vel,new2_points,Z,R_c2w,e);
   %% Threshold outputs in a range
   for r=1:length(Vel)
       if Vel(r)>2
           Vel(r)=2;
       elseif Vel(r)<-2
           Vel(r)=-2;
       end
  end
   %% Fix the linear velocity
   Omg=transpose(R_c2w)*Vel(4:6,1);
   vel=transpose(R_c2w)*(Vel(1:3,1)+cross(Vel(4:6,1),CTB));
   Vel=[vel;Omg];
    %% ADD SOME LOW PASS FILTER CODE
    filter=1;              % Change fil between 1 to 3
    if n>filter+1
        for r=1:filter-1
            Vel=Vel+estimated_V(:,n-r);
        end
    end
    Vel=Vel/filter;
    estimated_V(:,n)=Vel;
    Estimatedv(:,n)=Vel; 
   for r=1:length(Estimatedv(:,n))
       if Estimatedv(r,n)>2
           Estimatedv(r,n)=2;
       elseif Estimatedv(r,n)<-2
           Estimatedv(r,n)=-2;
       end
  end
    
end
    plotData(Estimatedv, sampledData, sampledVicon, sampledTime, datasetNum)
