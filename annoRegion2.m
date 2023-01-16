function Mask = annoRegion2(C1,C2,C3,C4,BW1,BW2)
% ANNOLINE Annotate image with lines
%
% Returns BW mask. Jumps randomly in z / switches ramdomly between xy xz yz
MP = get(0,'MonitorPositions'); % find position of second monitor to get figure window over there
Mask = BW1; % start Mask at BW1
close all
k=1;


%MP(2,:);

POS1 = [1           1        128         80];
POS2 = [1           -30        128         80];

while k>0
    if rand>0.666666 % annotate in xy plane
        
        rr = randi(size(C1,3)-2)+1;% select random slice
        
        sideString = 'XY';
        titleString = [sideString ' , zz = ' num2str(rr) '. Up and down is red and blue, green is present slice, purple is ZO1'];

        close all
      
        BW1rand = Mask(:,:,rr);
        BW2rand = BW2(:,:,rr);
        
        BWrand = BW1rand;
        BWrand(BW2rand==1)=2;
                
        COLVis = uint8(zeros(size(C1,1),size(C1,2),3)); 
        COLVis(:,:,1) = 2.*C1(:,:,rr-1)+ C3(:,:,rr);
        COLVis(:,:,2) = 2.*C1(:,:,rr) ;%+ C3(:,:,rr);
        COLVis(:,:,3) = 2.*C1(:,:,rr+1)+ C3(:,:,rr);
                
        COLVis2 = COLVis;% uint8(zeros(size(C1,1),size(C1,2),3)); 
                
       
        COLVis2(:,:,1) = max(4.*C4(:,:,rr-1:rr+1),[],3);
        COLVis2(:,:,2) = C2(:,:,rr); %+ 2.*C4(:,:,rr);
        COLVis2(:,:,3) = max(4.*C4(:,:,rr-1:rr+1),[],3);
        
        factor = 1;
        
        figure('OuterPosition',POS1)
        imshow(labeloverlay((factor.*COLVis),BWrand),'InitialMagnification',800); hold on; title(titleString)
        axis off
        set(gcf,'Visible','on')
        
        figure('OuterPosition',POS2)
        imshow(labeloverlay(factor.*COLVis2,BWrand),'InitialMagnification',800)
        axis off
        set(gcf,'Visible','on')
        
        lineMask = BW1rand;% zeros(size(BW1rand));%
        set(gca,'nextplot','replacechildren')
        
        h = drawpolygon(gca,'LineWidth',5);% a call to drawpolyline deletes previous polyline objects...
        % (double mouse press finishes line) and returns h
        lineMaskTemp = createMask(h);
        
        close all
        
        k = input('Continue drawing press 1, end press 0, undo press 2 :')
        
        if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
            lineMask(lineMaskTemp==1)=1;% ...  Add new line to it
        end
        
        %    imshow(labeloverlay(factor.*Srand,imdilate(lineMask,strel('disk',1))),'InitialMagnification',800)
        imshow(labeloverlay((factor.*squeeze(COLVis)),lineMask),'InitialMagnification',800); hold on; title('Up and down is red and blue, green is present slice, purple is ZO1')
        axis off
        
        Mask(:,:,rr) = lineMask;
        
    elseif rand>0.5 % annotate in zy plane

        rr = randi(size(C1,1)-2)+1;% select random slice

        sideString = 'ZY';
        titleString = [sideString ' , xx = ' num2str(rr) '. Up and down is red and blue, green is present slice, purple is ZO1'];

        
        close all
        
        BW1rand = permute(Mask(rr,:,:),[3 2 1]);
        BW2rand = permute(BW2(rr,:,:),[3 2 1]);
        
        BWrand = BW1rand;
        BWrand(BW2rand==1)=2;

                
        COLVis = uint8(zeros(size(C1,3),size(C1,2),3)); 
        COLVis(:,:,1) = 2.*permute(C1(rr-1,:,:),[3 2 1])+ permute(C3(rr,:,:),[3 2 1]);
        COLVis(:,:,2) = 2.*permute(C1(rr,:,:),[3 2 1]);% + permute(C3(rr,:,:),[3 2 1]);
        COLVis(:,:,3) = 2.*permute(C1(rr+1,:,:),[3 2 1])+ permute(C3(rr,:,:),[3 2 1]);
        
        COLVis2 = uint8(zeros(size(C1,3),size(C1,2),3)); 
%         COLVis2(:,:,1) = permute(C2(rr-1,:,:),[3 2 1])+ 2.*permute(C4(rr,:,:),[3 2 1]);
%         COLVis2(:,:,2) = permute(C2(rr,:,:),[3 2 1]) ;%+ 2.*permute(C4(rr,:,:),[3 2 1]);
%         COLVis2(:,:,3) = permute(C2(rr+1,:,:),[3 2 1])+ 2.*permute(C4(rr,:,:),[3 2 1]);

        COLVis2(:,:,1) = max( 4.*permute(C4(rr-1:rr+1,:,:),[3 2 1]),[],3);
        COLVis2(:,:,2) = permute(C2(rr,:,:),[3 2 1]) ;%+ 2.*permute(C4(rr,:,:),[3 2 1]);
        COLVis2(:,:,3) = max( 4.*permute(C4(rr-1:rr+1,:,:),[3 2 1]),[],3);
        
        factor = 1;
        
        
        figure('OuterPosition',[1           1        1280         800])
        imshow(labeloverlay(factor.*COLVis,BWrand),'InitialMagnification',800); hold on; title(titleString)
        axis off
        set(gcf,'Visible','on')
        
        figure('OuterPosition',POS2)
        imshow(labeloverlay(factor.*COLVis2,BWrand),'InitialMagnification',800)
        axis off
        set(gcf,'Visible','on')
        
        lineMask = BW1rand;%
        
        set(gca,'nextplot','replacechildren')
        
        h = drawpolygon(gca,'LineWidth',5);% a call to drawpolyline deletes previous polyline objects...
        % (double mouse press finishes line) and returns h
        lineMaskTemp = createMask(h);
        
        
        close all
        
        k = input('Continue drawing press 1, end press 0, undo press 2 :')
        
        if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
            lineMask(lineMaskTemp==1)=1;% ...  Add new line to it
        end
        
        axis off
        
        
        
        Mask(rr,:,:) = permute(lineMask,[3 2 1]);
                
        
        
        
    else % annotate in xz plane
        rr = randi(size(C1,2)-2)+1;% select random slice        
        sideString = 'ZX';
        titleString = [sideString ' , yy = ' num2str(rr) '. Up and down is red and blue, green is present slice , purple is ZO1'];

        close all
        
        BW1rand = permute(Mask(:,rr,:),[3 1 2]);
        BW2rand = permute(BW2(:,rr,:),[3 1 2]);
        
        BWrand = BW1rand;
        BWrand(BW2rand==1)=2;

                
        COLVis = uint8(zeros(size(C1,3),size(C1,2),3)); 
        COLVis(:,:,1) = 2.*permute(C1(:,rr-1,:),[3 1 2])+ permute(C3(:,rr,:),[3 1 2]);
        COLVis(:,:,2) = 2.*permute(C1(:,rr,:),[3 1 2]);% + permute(C3(:,rr,:),[3 1 2]);
        COLVis(:,:,3) = 2.*permute(C1(:,rr+1,:),[3 1 2])+ permute(C3(:,rr,:),[3 1 2]);

        
        COLVis2 = uint8(zeros(size(C1,3),size(C1,2),3)); 
%         COLVis2(:,:,1) = permute(C2(:,rr-1,:),[3 1 2])+ 2.*permute(C4(:,rr,:),[3 1 2]);
%         COLVis2(:,:,2) = permute(C2(:,rr,:),[3 1 2]);% + 2.*permute(C4(:,rr,:),[3 1 2]);
%         COLVis2(:,:,3) = permute(C2(:,rr+1,:),[3 1 2])+ 2.*permute(C4(:,rr,:),[3 1 2]);

        COLVis2(:,:,1) =  max(4.*permute(C4(:,rr-1:rr+1,:),[3 1 2]),[],3);
        COLVis2(:,:,2) = permute(C2(:,rr,:),[3 1 2]);% + 2.*permute(C4(:,rr,:),[3 1 2]);
        COLVis2(:,:,3) =  max(4.*permute(C4(:,rr-1:rr+1,:),[3 1 2]),[],3);
        
        factor = 1;

        figure('OuterPosition',POS1)
        imshow(labeloverlay((factor.*COLVis),BWrand),'InitialMagnification',800); hold on; title(titleString)
        axis off
        set(gcf,'Visible','on')
        
        figure('OuterPosition',POS2)
        imshow(labeloverlay(factor.*COLVis2,BWrand),'InitialMagnification',800)
        axis off
        set(gcf,'Visible','on')
        
        lineMask = BW1rand;
        
        set(gca,'nextplot','replacechildren')
        
        h = drawpolygon(gca,'LineWidth',5)% a call to drawpolyline deletes previous polyline objects...
        % (double mouse press finishes line) and returns h
        lineMaskTemp = createMask(h);
        
        
        close all
        
        k = input('Continue drawing press 1, end press 0, undo press 2 :')
        
        if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
            lineMask(lineMaskTemp==1)=1;% ...  Add new line to it
        end
        
        axis off
                
        Mask(:,rr,:) = permute(lineMask,[2 3 1]);
        
        
    end
    
end
end