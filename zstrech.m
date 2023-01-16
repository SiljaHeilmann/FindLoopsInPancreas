function [Iout] = zstrech(Iin)

[Xw, Yw, Zh] = size(Iin);


Iout = zeros(Xw,Yw,Zh*2-1);

for zz=2:Zh
    Iout(:,:,(zz-1)*2-1) = Iin(:,:,zz-1);   
    Iout(:,:,(zz-1)*2)   = mean(Iin(:,:,zz-1:zz),3);
    
end

end

