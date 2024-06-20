function De_CellCtr = HFBulb_DeCellCtr(cellctr, CellConnectivity, Cellid,  Ncell ) 

global Idx_Epi  Idx_Fan Idx_DP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
Mu_r = 15;
Mu_a = 3;
%Mu = 5;
% one cell wide
d0=0.5;


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
            
        
        %lind, a list of the indices of neighboring cells to ncell.
        for mcell = 1:1:length(lind)
            
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
            De_ncell = De_ncell + Mu*CellVec*(d0/norm(CellVec)-1);
        end
        
        De_CellCtr(ncell,:) = De_ncell ; 
    end

end

 

 
