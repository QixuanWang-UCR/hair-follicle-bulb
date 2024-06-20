function HFBulb_WT_MultiCol(nsim)
% clc; clf; close all; clear all;


%WT HF bulb
% Colonal drift of Mx
% Mx cells unlimited dividing potential

load StructHFB_Prep.mat

global Idx_Epi Idx_Fan Idx_DP  

Idx_Epi = 2;
Idx_DP = 1;
Idx_Fan = 0; 
x_midline = 4;

n_SaveData = length(StructHFB);

Ncell = StructHFB(n_SaveData).Ncell;
cellctr = StructHFB(n_SaveData).cellctr;
Cellid = StructHFB(n_SaveData).Cellid;
celltime = StructHFB(n_SaveData).celltime;
Time0 = StructHFB(n_SaveData).Time;

%Track cells' colony 
%originally colony 0
CellCol = zeros(1,Ncell);
count_col= 0;
N_Max_Col = 10;
N_Col_Time = 200;

%use yellow -> blue to show the colony number
%0 -> yellow [1,1,0]
%N_Max_Col -> cyan [0,1,1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NSaveTime = 500;
dt = 0.01;
NTime = N_Max_Col*N_Col_Time/dt;  



DivT_Avg = 50;
DivT_Std = 10;


Small_d= 0.1;

for ncell=1:1:Ncell
    if Cellid(ncell)==Idx_Epi  
       
        celltime(ncell)=DivT_Avg+DivT_Std*randn;

    end
end

yTop = 10;

%get the initial cell connectivity
TriMesh = delaunay(cellctr(:,1),cellctr(:,2));
TR = triangulation(TriMesh,cellctr(:,1),cellctr(:,2));
CellConnectivity = edges(TR);



 

for ntime = 1:1:NTime

    Time = ntime*dt + Time0;

%         StructHFB(n_SaveData).Ncell = Ncell;
%         StructHFB(n_SaveData).cellctr = cellctr;
%         StructHFB(n_SaveData).Cellid = Cellid;
%         StructHFB(n_SaveData).celltime = celltime;

    if mod(Time,N_Col_Time) == 1

         %Find the lowest pair of epi cells in the current system, mark
            %them as the next col
    
            count_col = count_col+1;

            cellctr_X = cellctr(:,1);
            cellctr_Y = cellctr(:,2);

            [~,min_left_y] = min(cellctr_Y((Cellid==Idx_Epi)&(cellctr_X<x_midline)));
            left_list = find((Cellid==Idx_Epi)&(cellctr_X<x_midline));
            min_left_id = left_list(min_left_y);
    
            CellCol(min_left_id) = count_col;
            celltime(min_left_id)=DivT_Avg+DivT_Std*randn;
          

            [~,min_right_y] = min(cellctr_Y((Cellid==Idx_Epi)&(cellctr_X>x_midline)));
            right_list = find((Cellid==Idx_Epi)&(cellctr_X>x_midline));
            min_right_id = right_list(min_right_y);
    
            CellCol(min_right_id) = count_col;
            celltime(min_right_id)=DivT_Avg+DivT_Std*randn;
    end
     
    %%%%%%%%% Cell Divide %%%%%%%%%%%%%%%%
  
    DeleteCell = [];
    CountDiv = 0;

    HighCell = [];
    CountHigh = 0;

    for ncell = 1:1:Ncell
        %Record the index of the dividing cells

        [row, col] = find(CellConnectivity==ncell);
        col = mod(col,2)+1;
    
        lind = sub2ind(size(CellConnectivity), row, col);
 
        AtchDP = 0;            
        
        %lind, a list of the indices of neighboring cells to ncell.
        for mcell = 1:1:length(lind)
            
            m_cell_id = CellConnectivity(lind(mcell));
                        
            %check if this neighbor is DP
            if Cellid(m_cell_id) == Idx_DP
                AtchDP = 1; 

            end

        end

 
        if Cellid(ncell)==Idx_Epi && celltime(ncell)<=0 && AtchDP == 1

            CountDiv = CountDiv + 1;

            %This cell is dividing, record it
            DeleteCell(CountDiv) = ncell;

            % Add daughter cells - position

            %if AtchDP == 1
                %attach to some DP
                divideangle = 0;
            %else
                %not attach to DP
                %divideangle=pi*rand;
            %end
            dividevec=[cos(divideangle);sin(divideangle)].';

            cellctr(Ncell+2*CountDiv-1,:)=cellctr(ncell,:)+Small_d*dividevec;
            cellctr(Ncell+2*CountDiv,:)=cellctr(ncell,:)-Small_d*dividevec;
              
            % Update daughter cell info
            Cellid(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=Idx_Epi;
            CellCol(Ncell+2*CountDiv-1:Ncell+2*CountDiv) = CellCol(ncell);
            celltime(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=[DivT_Avg+DivT_Std*randn,DivT_Avg+DivT_Std*randn];


        elseif cellctr(ncell,2)>yTop
                %Remove epi cells that are too high
                CountHigh = CountHigh + 1;
                HighCell(CountHigh) = ncell;

        end
          
    end

    DeleteCell = [DeleteCell, HighCell];

    cellctr(DeleteCell,:)=[]; 
    Cellid(DeleteCell)=[];
    celltime(DeleteCell)=[];
    CellCol(DeleteCell) = [];
    Ncell = length(Cellid);

 
    %%%%%%%%%% Update Dynamics %%%%%%%%%%%%%%%%

    TriMesh = delaunay(cellctr(:,1),cellctr(:,2));
    TR = triangulation(TriMesh,cellctr(:,1),cellctr(:,2));
    CellConnectivity = edges(TR);

    De_CellCtr = HFBulb_DeCellCtr_1(cellctr, CellConnectivity, Cellid, Ncell );
    cellctr = cellctr + De_CellCtr*dt;

    %%%%%%%%%% Boundary conditions %%%%%%%%%%%%%%

    for ncell = 1:1:Ncell

        if Cellid(ncell) == Idx_Epi
            %No leaking walls
            if cellctr(ncell,2) > 2
                %left
                if cellctr(ncell,1)<1.2
                    cellctr(ncell,1) = 1.2;
                elseif cellctr(ncell,1)>6.8
                    cellctr(ncell,1) = 6.8;
                end
            else
                %bottom
                cell_R = sqrt((cellctr(ncell,1)-4)^2+ (cellctr(ncell,2)-2)^2);
                if cell_R > 2.9
                    cellctr(ncell,:) = [4,2]+ (cellctr(ncell,:)-[4,2])*2.9/cell_R;

                end
            end

            %No entering DP region
            if cellctr(ncell,2) <3.6

                if cellctr(ncell,1)>3.4 && cellctr(ncell,1)<4
                    cellctr(ncell,1) = 3.3;
                elseif cellctr(ncell,1)<4.7 && cellctr(ncell,1)>4
                    cellctr(ncell,1) = 4.7;
                end

            end

        end

    end
 

    %%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%
%   
    if mod(ntime, NSaveTime) == 0 
% 
%         figure(1) 
%         clf
%         
%         
%         [Vor_v,Vor_c] = voronoin(cellctr);
%         
%         
%         for ncell = 1:1:Ncell
%         
%           if Cellid(ncell)==Idx_Epi   
%         
%               n_col = CellCol(ncell);
%               d_col = 1/N_Max_Col;
%         
%               patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2), [1,1,0] + n_col*[-d_col,0,d_col]);
%               hold on
%         
%           elseif Cellid(ncell)==Idx_DP 
%         
%               patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'k');
%               hold on
%         
%           end
%         
%         end
%         
%         plot(cellctr(Cellid==Idx_Fan,1),cellctr(Cellid==Idx_Fan,2),'b.','MarkerSize',10);
%         hold on
%         
%         Time = ntime*dt;
%         
%         StrTitle = sprintf('Time = %g', Time);
%         title(StrTitle)
%         
%         axis([-2,10,-2,16])
%         axis equal

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Save initial data
        % 
   
        n_SaveData = n_SaveData + 1;

        StructHFB(n_SaveData).Time = Time;
        StructHFB(n_SaveData).Ncell = Ncell;
        StructHFB(n_SaveData).cellctr = cellctr;
        StructHFB(n_SaveData).Cellid = Cellid;
        StructHFB(n_SaveData).celltime = celltime;
        StructHFB(n_SaveData).CellCol = CellCol;

    end
   
    celltime = celltime - dt;



    
end
 

StructTitleStr = ['BasDivZeroANgle_WT_MultiColT2000' num2str(nsim) '.mat']; 
save(StructTitleStr,'StructHFB','-v7.3')

 


       