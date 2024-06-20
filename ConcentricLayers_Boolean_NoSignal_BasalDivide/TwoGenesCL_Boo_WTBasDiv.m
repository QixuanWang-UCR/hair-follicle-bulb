%WT HF bulb
function TwoGenesCL_Boo_WTBasDiv(s_str, NBool)

% clf; clear all; clc

load StructHFB_Prep_Boo.mat

global Idx_Epi Idx_Fan Idx_DP  DP_b DP_t d0

Idx_Epi = 2;
Idx_DP = 1;
Idx_Fan = 0; 
% Idx_Blue = 3;

d0=0.5;

n_SaveData = length(StructHFB);

Ncell = StructHFB(n_SaveData).Ncell;
cellctr = StructHFB(n_SaveData).cellctr;
Cellid = StructHFB(n_SaveData).Cellid;
celltime = StructHFB(n_SaveData).celltime;

Time0 = StructHFB(n_SaveData).Time;
 CellGene =-1*ones(Ncell,2);
% 
% figure(1) 
% clf
% 
% [Vor_v,Vor_c] = voronoin(cellctr);
%  
% for ncell = 1:1:Ncell
% 
%     if Cellid(ncell)==Idx_Epi  
% 
%         patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'y');
%         hold on
% % 
% %     elseif Cellid(ncell)==Idx_Blue
% %         
% %         patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'c');
% %         hold on
% 
%     elseif Cellid(ncell)==Idx_DP 
% 
%         patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'k');
%         hold on
% 
%     end
% 
% end
% 
% 
% plot(cellctr(Cellid==Idx_Fan,1),cellctr(Cellid==Idx_Fan,2),'b.','MarkerSize',10);
% hold on
% plot(cellctr(Cellid==Idx_DP,1),cellctr(Cellid==Idx_DP,2),'k.','MarkerSize',10);
% hold on

% axis([-2,10,-2,16])
% axis equal

% StrTitle = sprintf('Time = %g', 0);
% title(StrTitle)

% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NFigTime = 1000;
NSaveData = 100;
NTime =200000;
dt = 0.01;


DivT_Avg = 50;
DivT_Std = 10;

Small_d= 0.1;




for ncell=1:1:Ncell
    if Cellid(ncell)==Idx_Epi 
       
        celltime(ncell)=DivT_Avg+DivT_Std*randn;

    end
end

yTop = 12;

DP_b = min(cellctr(Cellid==Idx_DP,2));
DP_t = max(cellctr(Cellid==Idx_DP,2))+d0;
 
TriMesh = delaunay(cellctr(:,1),cellctr(:,2));
TR = triangulation(TriMesh,cellctr(:,1),cellctr(:,2));
CellConnectivity = edges(TR);
% tic 

 

for ntime = 1:1:NTime

    Time = ntime*dt + Time0;

  


    %%%%%%%%% Cell Divide %%%%%%%%%%%%%%%%
  
    DeleteCell = [];
    CountDiv = 0;

    HighCell = [];
    CountHigh = 0;

    for ncell = 1:1:Ncell

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


        %Record the index of the dividing cells
 
        if Cellid(ncell)==Idx_Epi && celltime(ncell)<=0 && AtchDP == 1

            CountDiv = CountDiv + 1;

            %This cell is dividing, record it
            DeleteCell(CountDiv) = ncell;

            % Add daughter cells - position

            divideangle=0;
            dividevec=[cos(divideangle);sin(divideangle)].';

            cellctr(Ncell+2*CountDiv-1,:)=cellctr(ncell,:)+Small_d*dividevec;
            cellctr(Ncell+2*CountDiv,:)=cellctr(ncell,:)-Small_d*dividevec;
              
            % Update daughter cell info
            Cellid(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=Idx_Epi;
            celltime(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=[DivT_Avg+DivT_Std*randn,DivT_Avg+DivT_Std*randn];
            CellGene(Ncell+2*CountDiv-1,:)= CellGene(ncell,:);
            CellGene(Ncell+2*CountDiv,:)= CellGene(ncell,:);

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
    CellGene(DeleteCell,:)=[];
    Ncell = length(Cellid);

 
    %%%%%%%%%% Update Dynamics %%%%%%%%%%%%%%%%

    TriMesh = delaunay(cellctr(:,1),cellctr(:,2));
    TR = triangulation(TriMesh,cellctr(:,1),cellctr(:,2));
    CellConnectivity = edges(TR);

    [De_CellCtr, CellGene] = TwoGenes_DeCellCtr_3(cellctr, CellConnectivity, ntime, Cellid, s_str, Ncell, CellGene, NBool );
    cellctr = cellctr + De_CellCtr*dt;

    celltime = celltime - dt;

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
%     if mod(ntime, NFigTime) == 0 
%     %         every few step, draw
%     
%     
%             figure(1) 
%             clf
%     
%                 
%             [Vor_v,Vor_c] = voronoin(cellctr);
%      
%             for ncell = 1:1:Ncell 
%                  
%                       if Cellid(ncell)==Idx_Epi && max( CellGene(ncell,:) )< 0
%             
%                           patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'w');
%                           hold on
%             
%                       elseif Cellid(ncell)==Idx_Epi && max(CellGene(ncell,:) )>= 0
% 
%                           gHS = CellGene(ncell,1);
%                          
%                           gI = CellGene(ncell,2);
% 
%                           if gHS+gI > 1 || gHS+gI == 0
%                       
%                             patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'c');
%                             hold on
%                           elseif gHS == 1
%                               patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'r');
%                               hold on
%                           
%                           elseif gI == 1
%                               patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'g');
%                               hold on
%                           end
%             
%                       elseif Cellid(ncell)==Idx_DP 
%             
%                           patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'k');
%                           hold on
%             
%                       end
%             
%                   
%             
%                     plot(cellctr(Cellid==Idx_Fan,1),cellctr(Cellid==Idx_Fan,2),'b.','MarkerSize',10);
%                     hold on
%             
%                     Time = ntime*dt;
%             
%                     StrTitle = sprintf('Time = %g', Time);
%                     title(StrTitle)
%             
%                     axis equal
%                     axis([0,8,-2,10])
%             end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save data
    if mod(ntime, NSaveData) == 0
        n_SaveData = n_SaveData + 1;

        StructHFB(n_SaveData).Time = Time;
        StructHFB(n_SaveData).Ncell = Ncell;
        StructHFB(n_SaveData).cellctr = cellctr;
        StructHFB(n_SaveData).Cellid = Cellid;
        StructHFB(n_SaveData).celltime = celltime;
       StructHFB(n_SaveData).CellGene = CellGene;
    end
           
            
        
    
        


    
end

% toc

StructTitleStr = ['TwoGenesWT_ZeroAngle_Boo_NBool_' num2str(NBool) '_SStr_' num2str(s_str) '.mat'];
save(StructTitleStr,'StructHFB','-v7.3')

 


       