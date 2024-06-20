function [De_CellCtr, CellGene] = TwoGenes_DeCellCtr_3(cellctr, CellConnectivity, ntime, Cellid, s_str, Ncell, CellGene, NBool )


global Idx_Epi  Idx_Fan Idx_DP DP_b DP_t d0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
Mu_r = 15;
Mu_a = 3;
%Mu = 5;
% one cell wide
% s_str = 1;
% NBool = 50;

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
De_CellCtr = zeros(Ncell,2); 
Temp_CellGene = CellGene;
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
            
 %Count how many neighbor epi cells
        N_nb_epi = 0;

        %Count how many neighbor epi cells express the gene.
        pHS = 0;
        pI = 0;

        %lind, a list of the indices of neighboring cells to ncell.
        for mcell = 1:1:length(lind)
            
            m_cell_id = CellConnectivity(lind(mcell));
            CellVec = - cellctr(m_cell_id,:) + cellctr(ncell,:);

            %Check this neighbor's ID, Cellid(m_cell_id)
            %If DP -> Attach2DP = 1;
            %Else (do nothing)
            if Cellid(m_cell_id)==Idx_DP
                Attach2DP = 1;
             elseif Cellid(m_cell_id) == Idx_Epi && max(CellGene(m_cell_id,:))>0
                %If neighbor is an active Epi cell

                N_nb_epi = N_nb_epi + 1;

                if sum(CellGene(m_cell_id,:)) == 1
                    pHS = pHS + max(0, CellGene(m_cell_id,1));
                    
                    pI = pI + max(0, CellGene(m_cell_id,2));
                end
            end
            
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
        end
        
        De_CellCtr(ncell,:) = De_ncell ; 
        if Cellid(ncell)==Idx_Epi
            
            if Attach2DP == 1
                %Update gene Bool for epi cells attached to DP, based on their
                %y-axis
                cell_clr = min(1, max(0,(cellctr(ncell,2) - DP_b)/(DP_t - DP_b)));
    
                if cell_clr < 1/2 %IRS
    
                    Temp_CellGene(ncell,:) = [0,1];
    
               
    
                else %HS
    
                    Temp_CellGene(ncell,:) = [1,0];
    
                end

            elseif max(CellGene(ncell,:))>0 %If an epi cell is not attached to DP, update Bool model for this cell

                if N_nb_epi == 0
                    %no useful information from neighbors

                    pHS = 0;
                    pI = 0;

                else

                    pHS=s_str*pHS/N_nb_epi;
                    
                    pI=s_str*pI/N_nb_epi;

                    %Update SHS, SI
                    rHS = rand;
                    if rHS<pHS
                        sHS=1;
                    else
                        sHS=0;

                    end

                    

                    rI=rand;
                    if rI<pI
                        sI=1;
                    else
                        sI=0;
                    end

                    gHS = CellGene(ncell,1);
                    
                    gI = CellGene(ncell,2);

                    %Bool functions, async
                    Bool_Order = randperm(2);
                    for n_gene = 1:2
                        gene_id = Bool_Order(n_gene);

                        switch gene_id

                            case 1

                                Temp_CellGene(ncell,1) = (gHS && (~gI)) || sHS;
                                gHS = Temp_CellGene(ncell,1);
%                                 Async, need to use the updated info.,not
%                                 the old CellGene.(Temp_CellGene is the
%                                 updated one)(one line above)


                            case 2

                                Temp_CellGene(ncell,2) = (gI && (~gHS)) || sI;
                                gI = Temp_CellGene(ncell,2);

                            
                        end

                    end
                    
  
                end



            end
        end

        
 

    end

end
if mod(ntime, NBool) == 0  
    CellGene = Temp_CellGene;
end

 
