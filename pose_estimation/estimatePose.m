function [position, orientation, R_c2w] = estimatePose(data, t)
%% CHANGE THE NAME OF THE FUNCTION TO estimatePose
% Please note that the coordinates for each corner of each AprilTag are
% defined in the world frame, as per the information provided in the
% handout. Ideally a call to the function getCorner with ids of all the
% detected AprilTags should be made. This function should return the X and
% Y coordinate of each corner, or each corner and the centre, of all the
% detected AprilTags in the image. You can implement that anyway you want
% as long as the correct output is received. A call to that function
% should made from this function.
    %% Input Parameter Defination
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset
    
    %% Output Parameter Defination
    
    % position = translation vector representing the position of the
    % drone(body) in the world frame in the current time, in the order XYZ
    
    % orientation = euler angles representing the orientation of the
    % drone(body) in the world frame in the current time, in the order XYZ
    
    %R_c2w = Rotation which defines camera to world frame
    x1=zeros(length(data(t).id),1);x2=zeros(length(data(t).id),1);x3=zeros(length(data(t).id),1);x4=zeros(length(data(t).id),1);
    y1=zeros(length(data(t).id),1);y2=zeros(length(data(t).id),1);y3=zeros(length(data(t).id),1);y4=zeros(length(data(t).id),1);
    ID=data(t).id;
    P_1=data(t).p1; P1=(P_1)'; x1_dash= P1(:,1); y1_dash=P1(:,2);
    P_2=data(t).p2; P2=(P_2)'; x2_dash= P2(:,1); y2_dash=P2(:,2);
    P_3=data(t).p3; P3=(P_3)'; x3_dash= P3(:,1); y3_dash=P3(:,2);
    P_4=data(t).p4; P4=(P_4)'; x4_dash= P4(:,1); y4_dash=P4(:,2);
    X_Temp=[x1_dash,x2_dash,x3_dash,x4_dash];
    Y_temp=[y1_dash,y2_dash,y3_dash,y4_dash];
    for i=1:length(data(t).id)
        res=getCorner(ID(1,i));
        x1(i,1)=res(1,1);x2(i,1)=res(2,1);x3(i,1)=res(3,1);x4(i,1)=res(4,1); 
        y1(i,1)=res(5,1);y2(i,1)=res(6,1);y3(i,1)=res(7,1);y4(i,1)=res(8,1);
    end
    X=[x1,x2,x3,x4];Y=[y1,y2,y3,y4];
    Var_A=zeros(1,9);
    for i=1:length(data(t).id)
    for j=1:4
        xi=X(i,j);
        yi=Y(i,j);
        Xi=X_Temp(i,j);
        Yi=Y_temp(i,j);
        Row_1=[xi, yi, 1, 0, 0, 0, -(Xi)*xi, -(Xi)*yi, -Xi];
        Row_2=[0, 0, 0, xi, yi, 1, -(Yi)*xi, -(Yi)*yi, -Yi];
        Var_A=[Var_A;Row_1;Row_2];
    end
   end
   A=Var_A(2:end,:);
   [~,~,V_temp]=svd(A);
   H_temp=V_temp(:,9);
   H=[H_temp(1:3)';H_temp(4:6)';H_temp(7:9)'];
   H=H*sign(H(3,3));
   K=[311.0520, 0, 201.8724;0, 311.3885, 113.6210;0, 0, 1];
   K_invH_Mat=K\H;
   R_from_vect=[K_invH_Mat(:,1), K_invH_Mat(:,2), cross(K_invH_Mat(:,1),K_invH_Mat(:,2))];
   [U,~,V]=svd(R_from_vect);
   R33=det(U*(V)'); 
   R_dash=[1, 0, 0;0, 1, 0;0, 0, (R33)];
   R=U*R_dash*(V)';
   T=K_invH_Mat(:,3)/norm(K_invH_Mat(:,1));
   TMatW_to_C=[(R)',-(R).'*T;[0,0,0,1]];
   TMatC_to_B=[0.7071, -0.7071, 0.0000, -0.0400;-0.7071, -0.7071, 0.0000, 0;0.0000, 0.0000, -1.0000, -0.0300;0, 0, 0, 1.0000];
   TMatW_to_B=(TMatW_to_C)*(TMatC_to_B);
   orientation_dash=(rotm2eul(TMatW_to_B(1:3,1:3),'ZYX'))';
   orientation=orientation_dash(end:-1:1);
   position=TMatW_to_B(1:3,end);
   R_c2w=TMatW_to_C(1:3,1:3)';    
end