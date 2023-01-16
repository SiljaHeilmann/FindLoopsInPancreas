function Mask = annoLineS(Im,BW)
% ANNOLINE Annotate image with lines
%
% Returns BW mask. Jumps randomly in z
MP = get(0,'MonitorPositions'); % find position of second monitor to get figure window over there
Mask = BW; % start Mask at BW
close all
k=1;
while k>0
    if rand>0.5 % annotate in xy plane
        close all
        
        rr = randi(size(Im,3)-2)+1;% select random slice
        Srand = Im(:,:,rr,1);
        BWrand = Mask(:,:,rr);
        
        ImVis = Im;
        
        ImVis(:,:,rr,1) = Im(:,:,rr,1)+Im(:,:,rr,2).*0.5;
        
        factor = 1;
        
        figure('OuterPosition',[1           1        1280         800])
        imshow(imoverlay((factor.*ImVis(:,:,rr-1:rr+1,1)),BWrand,'green'),'InitialMagnification',400); hold on; title('Up and down is red and blue, green is present slice and Muc1')
        axis off
        set(gcf,'Visible','on')
        
        figure('OuterPosition',MP(2,:))
        imshow(imoverlay(factor.*Srand,BWrand,'green'),'InitialMagnification',400)
        axis off
        set(gcf,'Visible','on')
        
        lineMask = BWrand;% zeros(size(BWrand));%
        set(gca,'nextplot','replacechildren')
        
        h = drawpolyline(gca,'LineWidth',5)% a call to drawpolyline deletes previous polyline objects...
        % (double mouse press finishes line) and returns h
        lineMaskTemp = createMask(h);
        
        close all
        
        k = input('Continue drawing press 1, end press 0, undo press 2 :')
        
        if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
            lineMask(lineMaskTemp==1)=1;% ...  Add new line to it
        end
        
        %    imshow(imoverlay(factor.*Srand,imdilate(lineMask,strel('disk',1)),'green'),'InitialMagnification',400)
        imshow(imoverlay((factor.*squeeze(ImVis(:,:,rr-1:rr+1,1))),lineMask,'magenta'),'InitialMagnification',400); hold on; title('Up and down is red and blue, green is present slice and Muc1')
        axis off
        
        Mask(:,:,rr) = lineMask;
        
    elseif rand>0.5 % annotate in yz plane
        close all
        
        rr = randi(size(Im,2)-2)+1;% select random slice
        Srand = Im(:,rr,:,1);
        BWrand = Mask(:,rr,:);
        
        ImVis = Im;
        
        ImVis(:,rr,:,1) = Im(:,rr,:,1)+Im(:,rr,:,2).*0.5;
        
        factor = 1;
        
        
        figure('OuterPosition',[1           1        1280         800])
        imshow(imoverlay((factor.*permute(squeeze(ImVis(:,rr-1:rr+1,:,1)),[3  1  2])),permute(BWrand,[3 1 2]),'green'),'InitialMagnification',400); hold on; title('Up and down is red and blue, green is present slice and Muc1')
        axis off
        set(gcf,'Visible','on')
        
        figure('OuterPosition',MP(2,:))
        imshow(imoverlay(factor.*permute(Srand,[3 1  2]),permute(BWrand,[3 1 2]),'green'),'InitialMagnification',400)
        axis off
        set(gcf,'Visible','on')
        
        lineMask = BWrand;% zeros(size(BWrand));%
        lineMask = squeeze(lineMask);
        
        set(gca,'nextplot','replacechildren')
        
        h = drawpolyline(gca,'LineWidth',5)% a call to drawpolyline deletes previous polyline objects...
        % (double mouse press finishes line) and returns h
        lineMaskTemp = createMask(h);
        
        
        close all
        
        k = input('Continue drawing press 1, end press 0, undo press 2 :')
        
        if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
            lineMask(lineMaskTemp'==1)=1;% ...  Add new line to it
        end
        
        %    imshow(imoverlay(factor.*Srand,imdilate(lineMask,strel('disk',1)),'green'),'InitialMagnification',400)
        %imshow(imoverlay((factor.*ImVis(:,rr-1:rr+1,:,1)),lineMask,'magenta'),'InitialMagnification',400); hold on; title('Up and down is red and blue, green is present slice and Muc1')
        axis off
        
        
        Mask(:,rr,:) = squeeze(lineMask);
        
    else % annotate in xz plane
        close all
        
        rr = randi(size(Im,1)-2)+1;% select random slice
        Srand = Im(rr,:,:,1);
        BWrand = Mask(rr,:,:);
        
        ImVis = Im;
        
        ImVis(rr,:,:,1) = Im(rr,:,:,1)+Im(rr,:,:,2).*0.5;
        
        factor = 1;
        
        
        figure('OuterPosition',[1           1        1280         800])
        imshow(imoverlay((factor.*permute(squeeze(ImVis(rr-1:rr+1,:,:,1)),[3  2  1])),permute(BWrand,[3  2 1]),'green'),'InitialMagnification',400); hold on; title('Up and down is red and blue, green is present slice and Muc1')
        axis off
        set(gcf,'Visible','on')
        
        figure('OuterPosition',MP(2,:))
        imshow(imoverlay(factor.*permute(Srand,[3 2  1 ]),permute(BWrand,[3 2 1]),'green'),'InitialMagnification',400)
        axis off
        set(gcf,'Visible','on')
        
        lineMask = BWrand;% zeros(size(BWrand));%
        lineMask = squeeze(lineMask);
        
        set(gca,'nextplot','replacechildren')
        
        h = drawpolyline(gca,'LineWidth',5)% a call to drawpolyline deletes previous polyline objects...
        % (double mouse press finishes line) and returns h
        lineMaskTemp = createMask(h);
        
        hej=999
        size(lineMaskTemp)
        size(lineMask)
        
        close all
        
        k = input('Continue drawing press 1, end press 0, undo press 2 :')
        
        if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
            lineMask(lineMaskTemp'==1)=1;% ...  Add new line to it
        end
        
        %    imshow(imoverlay(factor.*Srand,imdilate(lineMask,strel('disk',1)),'green'),'InitialMagnification',400)
        %imshow(imoverlay((factor.*ImVis(:,rr-1:rr+1,:,1)),lineMask,'magenta'),'InitialMagnification',400); hold on; title('Up and down is red and blue, green is present slice and Muc1')
        axis off
        
        
        Mask(rr,:,:) = squeeze(lineMask);
        
        
    end
    
end
end