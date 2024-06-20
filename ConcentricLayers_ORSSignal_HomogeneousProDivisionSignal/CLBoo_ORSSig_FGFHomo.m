clf; clear all; clc

load StructHFB_Prep_Boo.mat

global Idx_Epi Idx_Fan Idx_DP  DP_b DP_t d0

Idx_Epi = 2;
Idx_DP = 1;
Idx_Fan = 0; 
% Idx_Blue = 3;
d0 = 0.5;
n_SaveData = 1;


Ncell = StructHFB(n_SaveData).Ncell;
cellctr = StructHFB(n_SaveData).cellctr;
Cellid = StructHFB(n_SaveData).Cellid;
celltime = StructHFB(n_SaveData).celltime;
Time0 = StructHFB(n_SaveData).Time;
% CellClr = StructHFB(n_SaveData).CellClr;
CellGene =-1*ones(Ncell,2);
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
% 
%     elseif Cellid(ncell)==Idx_Blue
%         
%         patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'c');
%         hold on
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
% 
% % axis([-2,10,-2,16])
% % axis equal
% 
% StrTitle = sprintf('Time = %g', 0);
% title(StrTitle)

% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NFigTime = 500;
NTime = 100000;
NSaveData = 100;
dt = 0.01;


DivT_Avg = 50;%100;
DivT_Std = 10;%20;


Small_d= 0.1;



for ncell=1:1:Ncell
    if Cellid(ncell)==Idx_Epi 
       
        celltime(ncell)=DivT_Avg+DivT_Std*randn;

    end
end

yTop = 12;
% yb=-0.5;
% yt=4;
DP_b = min(cellctr(Cellid==Idx_DP,2));
DP_t = max(cellctr(Cellid==Idx_DP,2))+d0;

TriMesh = delaunay(cellctr(:,1),cellctr(:,2));
TR = triangulation(TriMesh,cellctr(:,1),cellctr(:,2));
CellConnectivity = edges(TR);

% tic 

SignFgf = zeros(Ncell,1);
StructHFB(n_SaveData).SignFgf = SignFgf;

SignORS = zeros(Ncell,1);
StructHFB(n_SaveData).SignORS = SignORS; 

for ntime = 1:1:NTime

    Time = ntime*dt + Time0;


     
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
 
        Attach2DP = 0;            
        
        %lind, a list of the indices of neighboring cells to ncell.
        for mcell = 1:1:length(lind)
            
            m_cell_id = CellConnectivity(lind(mcell));
                        
            %check if this neighbor is DP
            if Cellid(m_cell_id) == Idx_DP
                Attach2DP = 1; 

            end

        end

        if Cellid(ncell)==Idx_Epi   && celltime(ncell)<=0

            CountDiv = CountDiv + 1;

            %This cell is dividing, record it
            DeleteCell(CountDiv) = ncell;

            % Add daughter cells - position

            if Attach2DP == 1
                %attach to some DP
                divideangle = 0;
            else
                %not attach to DP
                divideangle=pi*rand;
            end
            dividevec=[cos(divideangle);sin(divideangle)].';

            cellctr(Ncell+2*CountDiv-1,:)=cellctr(ncell,:)+Small_d*dividevec;
            cellctr(Ncell+2*CountDiv,:)=cellctr(ncell,:)-Small_d*dividevec;
              
            % Update daughter cell info
            Cellid(Ncell+2*CountDiv-1:Ncell+2*CountDiv)= Cellid(ncell);
            celltime(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=[DivT_Avg+DivT_Std*randn,DivT_Avg+DivT_Std*randn];
            SignFgf(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=SignFgf(ncell)/2;
            SignORS(Ncell+2*CountDiv-1:Ncell+2*CountDiv)=SignORS(ncell)/2;
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
    SignFgf(DeleteCell)=[];
    SignORS(DeleteCell)=[];
    Ncell = length(Cellid);
    CellGene(DeleteCell,:)=[];

 
    %%%%%%%%%% Update Dynamics %%%%%%%%%%%%%%%%

    TriMesh = delaunay(cellctr(:,1),cellctr(:,2));
    TR = triangulation(TriMesh,cellctr(:,1),cellctr(:,2));
    CellConnectivity = edges(TR);
    [Vor_v,Vor_c] = voronoin(cellctr);
    
    %Gene regulation by Boolean model
    [De_CellCtr, De_SignFgf, De_SignORS, CellGene] = CLBoo_DeCellCtr_ORSSig_FGFHomo(cellctr, CellConnectivity, Cellid, SignFgf,SignORS, Ncell,Vor_c, Vor_v, CellGene);

    cellctr = cellctr + De_CellCtr*dt;

    SignFgf = SignFgf + De_SignFgf*dt;
    SignFgf = max(0, SignFgf);
    SignFgf(cellctr(:,2)>11) = 0;

    SignORS = SignORS + De_SignORS*dt;
    SignORS = max(0, SignORS);
    SignORS(cellctr(:,2)>11) = 0;

    celltime = celltime - 0.5*(SignFgf.')*dt;

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
  if mod(ntime, NFigTime) == 0 
%     %         every few step, draw
%     
%     
            figure(1) 
        clf
    
        [Vor_v,Vor_c] = voronoin(cellctr);
    
        for ncell = 1:1:Ncell
    
            if Cellid(ncell)==Idx_Epi && max( CellGene(ncell,:) )< 0
    
              patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'w');
              hold on
    
           elseif Cellid(ncell)==Idx_Epi && max(CellGene(ncell,:) )>= 0

              gHS = CellGene(ncell,1);
             
              gI = CellGene(ncell,2);

              if gHS+gI > 1 || gHS+gI == 0
          
                patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'c');
                hold on
              elseif gHS == 1
                  patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'r');
                  hold on
              
              elseif gI == 1
                  patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'b');
                  hold on
              end
    
            elseif Cellid(ncell)==Idx_DP 
    
                patch(Vor_v(Vor_c{ncell},1),Vor_v(Vor_c{ncell},2),'k');
                hold on
    
            end
    
        end
    
        shp = alphaShape(cellctr(Cellid==Idx_Fan,1),cellctr(Cellid==Idx_Fan,2));
        plot(shp,'EdgeColor','w', 'FaceColor', 'w')

        %Time scale, with signal 0.5 day ~ T = 25
    
        TDay = (0.5/25)*Time;
         
         
    
        StrTitle = sprintf('%g day', 0.1*floor(10*TDay));
        title(StrTitle)
     
 
        axis equal
    
     
        axis([0,8,-1,9])
        set(gca,'FontSize',15);

        
        figure(2) 
        clf
        
        for ncell=1:1:Ncell
             if Cellid(ncell) == Idx_Fan
             else
                Cellclr=ColorMap_SignORS(SignORS(ncell));
                patch('XData',Vor_v(Vor_c{ncell},1),'YData',Vor_v(Vor_c{ncell},2),'EdgeColor',Cellclr,'FaceColor',Cellclr);
                hold on
             end
        end
        
        plot(shp,'EdgeColor','w', 'FaceColor', 'w')
        
        StrTitle = sprintf('%g day', floor(TDay));% 0.1*floor(10*TDay));
        title(StrTitle)
        
        
        axis equal
        
        
        axis([0,8,-1,9])
        set(gca,'FontSize',15);

 
 
    end
    if mod(ntime, NSaveData) == 0  
        n_SaveData = n_SaveData + 1;
 

        StructHFB(n_SaveData).Time = Time;
        StructHFB(n_SaveData).Ncell = Ncell;
        StructHFB(n_SaveData).cellctr = cellctr;
        StructHFB(n_SaveData).Cellid = Cellid;
        StructHFB(n_SaveData).celltime = celltime;
        StructHFB(n_SaveData).SignFgf = SignFgf;
        StructHFB(n_SaveData).SignORS = SignORS;
        StructHFB(n_SaveData).CellGene = CellGene;
        
     end
end
    
% end

% toc


StructTitleStr = ['CLBoo_ORSSog_FGFHomo.mat'];
save(StructTitleStr,'StructHFB','-v7.3')


function Cellclr=ColorMap_SignORS(c)
CM=jet(51);
Cmax=5;
d_c=Cmax/50;
Cscale=floor(min(Cmax,c)/d_c)+1;
Cellclr=CM(Cscale,:);
end
 
