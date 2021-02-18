function [Vel_Final] = velocityRANSAC(optV,optpos,Z,~,e)
%% CHANGE THE NAME OF THE FUNCTION TO velocityRANSAC
    %% Input Parameter Description
    % optV = The optical Flow
    % optPos = Position of the features in the camera frame 
    % Z = Height of the drone
    % R_c2w = Rotation defining camera to world frame
    % e = RANSAC hyper parameter    
    %% Output Parameter Description
    % Vel = Linear velocity and angualr velocity vector
    M=3;
    prob= 0.99;
    k=round(log10(1-prob)/log10(1-e.^M));
    A_Mat=[]; B_Vel=[];
     for i=1:length(optpos)
        A_row_1=[-1, 0, optpos(1,i)];
        A_row_2=[0, -1, optpos(2,i)];
        B_row_1=[optpos(1,i)*optpos(2,i), -(1+(optpos(1,i)).^2), optpos(2,i)];
        B_row_2=[1+(optpos(2,i)).^2, -optpos(1,i)*optpos(2,i), -optpos(1,i)];
        A_dash=[A_row_1/Z(:,i),A_row_2];
        B_dash=[B_row_1/Z(:,i),B_row_2];
        A_Mat=[A_Mat; A_dash; B_dash];
        B_Vel=[B_Vel; optV(:,i)];
     end
     index_dash=0; 
    for M=1:k       
        pr=randperm(length(optpos),3);
        C_Mat=[]; D_Vel=[]; 
        for j=1:length(pr)          
            C_row_1=[-1, 0, optpos(1,pr(j))];
            C_row_2=[0, -1, optpos(2,pr(j))];
            D_row_1=[optpos(1,pr(j))*optpos(2,pr(j)), -(1+(optpos(1,pr(j))).^2), optpos(2,pr(j))];
            D_row_2=[1+(optpos(2,pr(j))).^2, -optpos(1,pr(j))*optpos(2,pr(j)), -optpos(1,pr(j))];
            C_dash=[C_row_1/Z(:,pr(j)),C_row_2];
            D_dash=[(D_row_1)/Z(:,pr(j)),D_row_2];
            C_Mat=[C_Mat; C_dash; D_dash];
            D_Vel=[D_Vel; optV(:,pr(j))];
        end
        R_V=inv(C_Mat)*D_Vel;
        ER=A_Mat*R_V-B_Vel;
        index=0;       
        ID=[];
        for q=1:length(optV)
            error=sqrt((ER(2*q))^2+(ER((2*q)-1))^2);
            if error<0.015                      % Change error limit between 0.001 to 0.8
                index=index+1;
                ID=[ID;q];
            end
        end
        if index>index_dash
            index_dash=index;
            ID_Final=ID;
        end
    end
     E_Mat=[]; F_Vel=[];
        for len=1:length(ID_Final)
        E_row_1=[-1, 0, optpos(1,ID_Final(len))];
        E_row_2=[0, -1, optpos(2,ID_Final(len))];
        F_row_1=[optpos(1,ID_Final(len))*optpos(2,ID_Final(len)), -(1+(optpos(1,ID_Final(len))).^2), optpos(2,ID_Final(len))];
        F_row_2=[1+(optpos(2,ID_Final(len))).^2, -optpos(1,ID_Final(len))*optpos(2,ID_Final(len)), -optpos(1,ID_Final(len))];
        E_dash=[E_row_1/Z(:,ID_Final(len)),E_row_2];
        F_dash=[F_row_1/Z(:,ID_Final(len)),F_row_2];
        E_Mat=[E_Mat; E_dash; F_dash];
        F_Vel=[F_Vel; optV(:,ID_Final(len))];             
        end
           Vel_Final=E_Mat\F_Vel;
end