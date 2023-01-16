function Mask = addRegion(C1,C2,BW,angle,start_rr)
% ANNOLINE Annotate image with lines
%
% Returns BW mask. Jumps randomly in z / switches ramdomly between xy xz yz
MP = get(0,'MonitorPositions'); % find position of second monitor to get figure window over there
Mask = BW; % start Mask at BW
close all
k=1;

alpha = 0.8;



POS1 = [1           1        128         80];
POS2 = [1           -30      128         80];

while k>0
    switch angle
        case {'XY','xy','YX','yx'}
            % if rand>1;%0.3333333 % annotate in xy plane
            rr = start_rr;%size(C1,3);
            
            for ii = 1:size(C1,3)
                if rr<=size(C1,3)
                    %       rr = randi(size(C1,3)-2)+1;% select random slice
                    
                    sideString = 'XY';
                    titleString = [sideString ' , zz = ' num2str(rr) ''];
                    
                    close all
                    
                    L = bwlabeln(Mask);
                    
                    COLVis = labeloverlay3(C2.*2,L,alpha);% uint8(zeros(size(C1,1),size(C1,2),3));
                    
                    MaskOutline = bwperim(Mask);
                    MaskOutline2 = max(MaskOutline(:,:,rr-1:rr+1),[],3);
           
                    IM = COLVis(:,:,:,rr);
                    IM(:,:,1)=COLVis(:,:,1,rr)+uint8(MaskOutline2.*80)+C2(:,:,rr+1);
                    IM(:,:,2)=COLVis(:,:,2,rr)+uint8(MaskOutline2.*80)+C2(:,:,rr-1);
                    IM(:,:,3)=COLVis(:,:,3,rr)+uint8(MaskOutline2.*80)+8.*C1(:,:,rr);
                    
                    figure('OuterPosition',POS1)
                    imshow(IM,'InitialMagnification',800); hold on; title(titleString)
                    axis off
                    set(gcf,'Visible','on')
                    
                    %         figure('OuterPosition',POS2)
                    %         imshow(C2(:,:,rr),'InitialMagnification',800)
                    %         axis off
                    %         set(gcf,'Visible','on')
                    
                    set(gca,'nextplot','replacechildren')
                    
                    h = drawpolygon(gca,'LineWidth',5);% a call to drawpolyline deletes previous polyline objects...
                    % (double mouse press finishes line) and returns h
                    lineMaskTemp = createMask(h);
                    
                    close all
                    
                    k = input('Continue drawing press 1, end press 0, undo press 2, stay in slice press 3 :')
                    
                    if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
                        TempMask = zeros(size(Mask));
                        TempMask(:,:,rr) = lineMaskTemp;
                        
         %      TempMask = imdilate(TempMask,strel('sphere',1));
                        Mask(TempMask==1)=1;
                        
                    elseif k==2
                        k=0
                        break
                        
                        
                    elseif k==3
                        %  lineMask(imdilate(lineMaskTemp==1,strel('disk',1))==1)=0;% ...  remove
                        TempMask = zeros(size(Mask));
                        TempMask(:,:,rr) = lineMaskTemp;
                        
         %      TempMask = imdilate(TempMask,strel('sphere',1));
                        Mask(TempMask==1)=1;
                        rr=rr-1
                    end
                    
                    rr=rr+1;
                    
                    
                    %    elseif rand>1 % annotate in zy plane
                end
            end
            
        case {'ZY','YZ','zy','yz'}
            rr = start_rr;
            
            for ii = 1:size(C1,1)
                if rr<=size(C1,1)
                    
                    %        rr = randi(size(C1,1)-2)+1;% select random slice
                    
                    sideString = 'ZY';
                    titleString = [sideString ' , xx = ' num2str(rr) ''];
                    
                    
                    close all
                    
                    Maskrot = permute(Mask,[3 2 1]);
                    C1rot  = permute(C1,[3 2 1]);
                    C2rot = permute(C2,[3 2 1]);
                    
                    L = bwlabeln(Maskrot);
                    
                    COLVis = labeloverlay3(C2rot.*2,L,alpha);%permute(labeloverlay3(C2rot.*4,L,alpha),[4 2 3 1]);
                    
                    MaskOutline = bwperim(Maskrot);
                    MaskOutline2 = max(MaskOutline(:,:,rr-1:rr+1),[],3);
           
                    IM = COLVis(:,:,:,rr);
                    IM(:,:,1)=COLVis(:,:,1,rr)+uint8(MaskOutline2.*80)+C2rot(:,:,rr+1);
                    IM(:,:,2)=COLVis(:,:,2,rr)+uint8(MaskOutline2.*80)+C2rot(:,:,rr-1);
                    IM(:,:,3)=COLVis(:,:,3,rr)+uint8(MaskOutline2.*80)+8.*C1rot(:,:,rr);
                    
                    
                    
                    
                    figure('OuterPosition',[1           1        1280         800])
                    imshow(IM,'InitialMagnification',800); hold on; title(titleString)
                    axis off
                    set(gcf,'Visible','on')
                    
                    
                    lineMask = Maskrot;%
                    
                    set(gca,'nextplot','replacechildren')
                    
                    h = drawpolygon(gca,'LineWidth',5);% a call to drawpolyline deletes previous polyline objects...
                    % (double mouse press finishes line) and returns h
                    lineMaskTemp = createMask(h);
                    
                    
                    close all
                    
                    k = input('Continue drawing press 1, end press 0, undo press 2, stay in slice press 3 :')
                    
                    if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
                        % lineMask(imdilate(lineMaskTemp==1,strel('disk',1))==1)=0;% ...  Add new line to it
                        
                        TempMask = zeros(size(Mask));
                        TempMask(rr,:,:) = permute(imdilate(lineMaskTemp==1,strel('disk',1)),[3 2 1]);
                        
         %      TempMask = imdilate(TempMask,strel('sphere',1));
                        Mask(TempMask==1)=1;
                        
                    elseif k==2
                        k=0
                        break
                        
                    elseif k==3
                        
                        TempMask = zeros(size(Mask));
                        TempMask(rr,:,:) = permute(imdilate(lineMaskTemp==1,strel('disk',1)),[3 2 1]);
                        
         %      TempMask = imdilate(TempMask,strel('sphere',1));
                        Mask(TempMask==1)=1;
                        
                        rr=rr-1
                    end
                    
                    axis off
                    rr=rr+1;
                    
                end
            end
        case {'ZX','XZ','zx','xz'}
            
            rr = start_rr;
            for ii = 1:size(C1,2)
                if rr<=size(C1,2)
                    
                    
                    %rr = randi(size(C1,2)-2)+1;% select random slice
                    sideString = 'ZX';
                    titleString = [sideString ' , yy = ' num2str(rr) '. Up and down is red and blue, green is present slice , purple is ZO1'];
                    
                    close all
                    
                    
                    L = bwlabeln(Mask);
                    clear COLVis
                    COLVis = permute(labeloverlay3(C2.*4,L,alpha),[4 1 3 2]);
                    
                    
                    figure('OuterPosition',POS1)
                    imshow(COLVis(:,:,:,rr),'InitialMagnification',800); hold on; title(titleString)
                    axis off
                    set(gcf,'Visible','on')
                    
                    
                    set(gca,'nextplot','replacechildren')
                    
                    h = drawpolygon(gca,'LineWidth',5)% a call to drawpolyline deletes previous polyline objects...
                    % (double mouse press finishes line) and returns h
                    lineMaskTemp = createMask(h);
                    
                    
                    close all
                    
                    k = input('Continue drawing press 1, end press 0, undo press 2, stay in slice press 3 :')
                    
                    if k==0 || k==1 % if you press 2 you stay in loop but last line is not added
                        % lineMask(imdilate(lineMaskTemp==1,strel('disk',1))==1)=0;% ...  Add new line to it
                        
                        TempMask = zeros(size(Mask));
                        TempMask(:,rr,:) = permute(imdilate(lineMaskTemp==1,strel('disk',1)),[2 3 1]);
                        
         %      TempMask = imdilate(TempMask,strel('sphere',1));
                        Mask(TempMask==1)=1;
                    elseif k==2
                        k=0
                        break
                        
                        
                    elseif k==3
                        % lineMask(imdilate(lineMaskTemp==1,strel('disk',1))==1)=0;% ...  Add new line to it
                        
                        TempMask = zeros(size(Mask));
                        TempMask(:,rr,:) = permute(imdilate(lineMaskTemp==1,strel('disk',1)),[2 3 1]);
                        
         %      TempMask = imdilate(TempMask,strel('sphere',1));
                        Mask(TempMask==1)=1;
                        
                        
                        rr=rr-1
                    end
                    axis off
                    
                    rr=rr+1;
                end
            end
            
    end
end

hej =999

end
