function [De_CellCtr, De_SignFgf, CellClr] = HFBulb_DeDyn(cellctr, CellConnectivity, Cellid, SignFgf, Ncell,Vor_c, Vor_v, CellClr )

global Idx_Epi  Idx_Fan Idx_DP DP_b DP_t d0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
Mu_r = 15;
Mu_a = 3;
%Mu = 5;
% one cell wide
d0=0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
De_CellCtr = zeros(Ncell,2); 
De_SignFgf = zeros(Ncell,1);
alpha = 4;
deg = 0.5;

% Attach2DP = 0;
%update cell center positions
for ncell = 1:1:Ncell
    
    %For the current #ncell, find all neighboring cells index

    if  Cellid(ncell)==Idx_Fan 

    else %top cells
    
        [row, col] = find(CellConnectivity==ncell);
        col = mod(col,2)+1;
    
        lind = sub2ind(size(CellConnectivity), row, col);
 
        De_ncell = [0,0];
        De_sign_Fgf = 0;
        Lapf = 0;

        
        %lind, a list of the indices of neighboring cells to ncell.
        for mcell = 1:1:length(lind)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %cell movement
            
            m_cell_id = CellConnectivity(lind(mcell));
            CellVec = - cellctr(m_cell_id,:) + cellctr(ncell,:);


            
            %spring between neighboring cell centers
            if norm(CellVec) > d0
                %weak attachment
                Mu = Mu_a;
            else
                %strong repel
                Mu = Mu_r;
            end

            if Cellid(ncell)==Idx_DP
            else

            De_ncell = De_ncell + Mu*CellVec*(d0/norm(CellVec)-1);
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %signal


            if Cellid(m_cell_id)==Idx_Fan || cellctr(ncell,2)>11
            else
                commonVT=intersect(Vor_c{ncell},Vor_c{m_cell_id});
                if length(commonVT)==2
                    l=norm(Vor_v(commonVT(1),:)-Vor_v(commonVT(2),:));
                    Lapf = Lapf+(SignFgf(m_cell_id)-SignFgf(ncell))*l/norm(CellVec);

                end
            end

        end

        if cellctr(ncell,2)<=11

            vx=Vor_v(Vor_c{ncell},1);
            vy=Vor_v(Vor_c{ncell},2);
            Area=polyarea(vx,vy);
    
            if Cellid(ncell)==Idx_DP
                De_sign_Fgf = Lapf/Area + alpha*(1 - (cellctr(ncell,2) - DP_b)/(DP_t - DP_b));
            else
                De_sign_Fgf = Lapf/Area - deg*SignFgf(ncell) ; 
            end  

        end
        
        De_CellCtr(ncell,:) = De_ncell ; 
        De_SignFgf(ncell) = De_sign_Fgf;


        
    end

end

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
        if Attach2DP == 1
            CellClr(ncell) = min(1, max(0,(cellctr(ncell,2) - DP_b)/(DP_t - DP_b)));
        
        end
 end

 
