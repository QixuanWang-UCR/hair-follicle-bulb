function [De_CellCtr,De_SignORS, CellGene] = CLBoo_ORSSig_DeCellCtr(SignORS,cellctr, CellConnectivity,  Cellid,  Ncell, CellGene, Vor_c, Vor_v)


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

De_SignORS = zeros(Ncell,1);
alpha_ORS = 1;
deg_ORS = 1;
SigThrd = 0.1;


for ncell = 1:1:Ncell
    
    %For the current #ncell, find all neighboring cells index

    if  Cellid(ncell)==Idx_Fan 
    else %top cells
    
        [row, col] = find(CellConnectivity==ncell);
        col = mod(col,2)+1;
    
        lind = sub2ind(size(CellConnectivity), row, col);
 
        De_ncell = [0,0];

        De_sign_ORS = 0;
        Lapf = 0;

        Attach2DP = 0;
        Attach2Fan = 0;
            

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

             if Cellid(ncell)==Idx_DP
            else
                De_ncell = De_ncell + Mu*CellVec*(d0/norm(CellVec)-1);
             end

            if cellctr(ncell,2)>11
            elseif Cellid(m_cell_id)==Idx_Fan 
                Attach2Fan = 1;
            else
                commonVT=intersect(Vor_c{ncell},Vor_c{m_cell_id});
                if length(commonVT)==2
                    l=norm(Vor_v(commonVT(1),:)-Vor_v(commonVT(2),:));
                    Lapf = Lapf+(SignORS(m_cell_id)-SignORS(ncell))*l/norm(CellVec);
    
                end
            end
        end

        if cellctr(ncell,2)<=11

            vx=Vor_v(Vor_c{ncell},1);
            vy=Vor_v(Vor_c{ncell},2);
            Area=polyarea(vx,vy);
    

            De_sign_ORS = Lapf/Area - deg_ORS*SignORS(ncell) + alpha_ORS*Attach2Fan ; 
  

        end
        
        De_CellCtr(ncell,:) = De_ncell;
        De_SignORS(ncell) = De_sign_ORS;


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

          
                    pI= max(SignORS(ncell),SigThrd) -  SigThrd;

                    %Update SI

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

                                Temp_CellGene(ncell,1) = (gHS && (~gI)) ;
                                gHS = Temp_CellGene(ncell,1);


                            case 2

                                Temp_CellGene(ncell,2) = (gI && (~gHS)) || sI;
                                gI = Temp_CellGene(ncell,2);

                            
                        end

                    end
                    
  
                end



            end
        end

        
 
 
        %Check if connect to any of DP cell -> CellConnectivity
         
         %If attach to DP (If Attach2DP == 1), define the color based on its y:
         % CellClr(ncell) = (y-B)/(T-B)

         % %IF not, do nothing

%         if Attach2DP == 1
%             CellClr(ncell) = min(1, max(0,(cellctr(ncell,2) - DP_b)/(DP_t - DP_b)));
%         
%         end


end
 
CellGene = Temp_CellGene;
 

 
