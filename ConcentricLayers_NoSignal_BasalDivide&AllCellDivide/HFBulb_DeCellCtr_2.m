function [De_CellCtr, CellClr] = HFBulb_DeCellCtr_2(cellctr, CellConnectivity, Cellid, CellClr, Ncell )


global Idx_Epi  Idx_Fan Idx_DP DP_b DP_t d0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
Mu_r = 15;
Mu_a = 3;
%Mu = 5;
% one cell wide


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
De_CellCtr = zeros(Ncell,2); 

%update cell center positions
for ncell = 1:1:Ncell
    
    %For the current #ncell, find all neighboring cells index

    if  Cellid(ncell)==Idx_Fan || Cellid(ncell)==Idx_DP
    else %top cells
    
        [row, col] = find(CellConnectivity==ncell);
        col = mod(col,2)+1;
    
        lind = sub2ind(size(CellConnectivity), row, col);
 
        De_ncell = [0,0];

        Attach2DP = 0;
            


        %lind, a list of the indices of neighboring cells to ncell.
        for mcell = 1:1:length(lind)
            
            m_cell_id = CellConnectivity(lind(mcell));
            CellVec = - cellctr(m_cell_id,:) + cellctr(ncell,:);

            %Check this neighbor's ID, Cellid(m_cell_id)
            %If DP -> Attach2DP = 1;
            %Else (do nothing)
            if Cellid(m_cell_id)==Idx_DP
                Attach2DP = 1;
            
            end
            
            %spring between neighboring cell centers
            if norm(CellVec) > d0
                %weak attachment
                Mu = Mu_a;
            else
                %strong repel
                Mu = Mu_r;
            end
            De_ncell = De_ncell + Mu*CellVec*(d0/norm(CellVec)-1);
        end
        
        De_CellCtr(ncell,:) = De_ncell ; 

 
 
        %Check if connect to any of DP cell -> CellConnectivity
         
         %If attach to DP (If Attach2DP == 1), define the color based on its y:
         % CellClr(ncell) = (y-B)/(T-B)

         % %IF not, do nothing

        if Attach2DP == 1
            CellClr(ncell) = min(1, max(0,(cellctr(ncell,2) - DP_b)/(DP_t - DP_b)));
        
        end

    end

end

 

 
