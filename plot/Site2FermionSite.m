function FermionSiteNum=Site2FermionSite(GlobalSiteNum, Ly, Np)
% For SSH-Hubbard model
% Parameter: 
% - GlobalSiteNum: count all the sites, including fermion and boson, count
% from 1
% - Ly: vertical length, only count fermion sites
% - Np: The number of pesudosites
% Return: 
% - FermionSiterNum: only count fermion sites, count from 0

        
        %Ly>2 cases
        if(Ly>2)
            residue = mod(GlobalSiteNum, (2*Np+1)*Ly);
            if(residue<(Np+1)*Ly && mod(residue,Np+1)==0)
                FermionSiteNum=(GlobalSiteNum-residue)/(2*Np+1) + residue/(Np+1);
            else
                FermionSiteNum=-1;% Boson site case
            end
        elseif(Ly==2)
            residue = mod(GlobalSiteNum, 3*Np+2);
            if(residue==0)
                FermionSiteNum=GlobalSiteNum/(3*Np+2)*2;
            elseif(residue==Np+1)
                FermionSiteNum=GlobalSiteNum/(3*Np+2)*2+1;
            else
                FermionSiteNum=-1;
            end
        end
end