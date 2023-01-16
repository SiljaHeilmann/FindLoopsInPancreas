function [Lnew] = relabelForVolView(L)
% takes a labelled image and returns a labelled image where no index is
% above 127 (it relables regions with index between 1 and 127)

Lnew = uint16(L);
MAX = max(Lnew(:));
Lnew(Lnew>0)=0;


for ii = 1:max(1,ceil(MAX/127))
    TEMP1 = L>(ii-1)*127 & L<(ii)*127;% mask with all current regions having index between (ii-1)*127 and (ii)*127
    TEMP2 = uint16(L);
    TEMP2(TEMP1==0)=0;% isolate current regions 
    TEMP2 = TEMP2-(ii-1)*127;% Lower their indexes by subtracting 127. (note its uint16 so anything below zero will be zeros...).
    Lnew = Lnew+TEMP2;% add them in to Lnew
end


end