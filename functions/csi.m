function csi(filename, experimentName)
    if isempty(filename)
        quit(2);
    end

    data = load(filename);
    NeuP = data.NeuP;
    HeadX = data.HeadX;
    HeadY = data.HeadY;
    Speed = data.Speed;


    % Determine how many segments
    [seg,frame]=size(NeuP);
    seg = 1:seg;
    num_seg = numel(seg);
    
    % Change X dimension to an even number
    Xmin=fix(min(HeadX));
    Xmax=ceil(max(HeadX));
    if mod((Xmax-Xmin),2)>0
        Xmax=Xmax+1;
    end
    
    % Change Y dimension to an odd number
    Ymin=fix(min(HeadY));
    Ymax=ceil(max(HeadY));
    if mod((Ymax-Ymin),2)>0
        Ymax=Ymax+1;
    end
    
    % Default value: 5 cm/s
    minSpeed = 5;
    
    framesSlow=find(Speed<minSpeed);
    frames=find(Speed>=minSpeed);
    
    Pixelsize=2;
    
    PixelNoX=(Xmax-Xmin)/Pixelsize;
    PixelNoY=(Ymax-Ymin)/Pixelsize;
    
    z=length(HeadX);
    head=zeros(z,2);
    for k=1:z
        head(k,1)=HeadX(k);
        head(k,2)=HeadY(k);
    end

    SZOccu(PixelNoX,PixelNoY)=0;
    [SZ(:,1),~]=discretize(head(:,1),PixelNoX);
    [SZ(:,2),~]=discretize(head(:,2),PixelNoY);
    for m=1:PixelNoX
        for n=1:PixelNoY
            indices = find(SZ(frames,1)==m & SZ(frames,2)==n);
            [xx,~] = size(indices);
            SZOccu(m,n)=xx;
        end
    end

    meanP_SZOccu(PixelNoX,PixelNoY,numel(num_seg))=0;
    for i=1:num_seg
        for m=1:PixelNoX
            for n=1:PixelNoY
                indices = find(SZ(frames,1)==m & SZ(frames,2)==n);
                meanP_SZOccu(m,n,i)=(sum(NeuP(seg(i),indices))/numel(indices))*5;
            end
        end
    end

    meanP(numel(num_seg))=0;
    PSZOccu=SZOccu/(frame);

    for i=1:num_seg
        meanP(i)=(sum(NeuP(seg(i),frames))/numel(frames))*5;
    end

    MSI(PixelNoX,PixelNoY,numel(num_seg))=0;
    for i=1:num_seg
        for r=1:PixelNoX
            for s=1:PixelNoY
                MSI(r,s,i)=PSZOccu(r,s)*(meanP_SZOccu(r,s,i))*(log2((meanP_SZOccu(r,s,i)/meanP(i))));
                SPARtop(r,s,i)=(PSZOccu(r,s)*meanP_SZOccu(r,s,i));
                SPARbot(r,s,i)=(PSZOccu(r,s)*(meanP_SZOccu(r,s,i)^2));
            end
        end

        MSIsum = nansum(MSI(:,:,i),'all');
        MSIseg(i,1)=MSIsum;
        aa=nansum(SPARtop(:,:,i),'all');
        bb=nansum(SPARbot(:,:,i),'all');
        SPARsum = (aa^2)/bb;
        SPARseg(i,1)=SPARsum;
    end
    save(filename, "minSpeed", "framesSlow", "PixelNoX", "Pixelsize", "PixelNoY", "-append");
    %% Calculating CSI of NeuP for different behavioral vectors
    % similarity index: dot sum of two vectors divided by the product of their norm values.
    [~,frame]=size(NeuP);
    seg=1:numel(seg);
    num_seg = numel(seg);
    orderingSize = 5000;

    order = 1:8;
    allorders = perms(order);
    [items,~] = size(allorders);
    randorders = randperm(items,orderingSize);
    note_interval = fix(frame/8);
    notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Speed)];
    
    %shuffle behavioral vectors into 8 pieces
    order = 1:8;
    allorders = perms(order);
    [items,~] = size(allorders);
    randorders = randperm(items,orderingSize);
    note_interval = fix(frame/8);
    notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Speed)];
    
    simI_Speed(num_seg)=0;
    for i=1:num_seg
        simI_Speed(i)=dot(NeuP(seg(i),:),Speed)/(norm(NeuP(seg(i),:)*norm(Speed)));
    end

    % variables for shuffled behavioral vectors
    shufSpeed(orderingSize,frame)=0;
    % creating shuffled behavioral vectors
    for o=1:8
        tempSpeed.piece(o).seq=Speed((notes(o)+1):notes(o+1));
    end
    for k=1:orderingSize
        shufSpeed(k,:)=[[tempSpeed.piece(allorders(randorders(k),1)).seq];[tempSpeed.piece(allorders(randorders(k),2)).seq];[tempSpeed.piece(allorders(randorders(k),3)).seq];[tempSpeed.piece(allorders(randorders(k),4)).seq];[tempSpeed.piece(allorders(randorders(k),5)).seq];[tempSpeed.piece(allorders(randorders(k),6)).seq];[tempSpeed.piece(allorders(randorders(k),7)).seq];[tempSpeed.piece(allorders(randorders(k),8)).seq]];
    end
    
    % similarity index of randomized data: dot sum of two vectors divided by the product of their
    % norm values (5000 randomized behavioral vectors).
    simI_shufSpeed(num_seg,orderingSize)=0;
    for i=1:num_seg
        for j=1:orderingSize
            simI_shufSpeed(i,j)=dot(NeuP(seg(i),:),shufSpeed(j,:))/(norm(NeuP(seg(i),:)*norm(shufSpeed(j,:))));
        end
    end
    
    % Counting how many randomized data are smaller than actual data
    SSeg_Speed(num_seg)=0;
    for m=1:num_seg
        SSeg_Speed(m)=numel(find(simI_Speed(m)>simI_shufSpeed(m,:)));
    end
    suffixes = [""];
    if experimentName == "SDT"
        suffixes = ["_A", "_N"];
    end
    listCSI = struct();
    places = struct();
    for k=1:length(suffixes)
        suffix = suffixes(k);
        CSI = struct();
        place = struct();
        mouseName = "mouse" + suffix;
        DistanceH = data.(mouseName).DistanceH;
        Compass = data.(mouseName).Compass;
        Behav_A50_D10 = data.(mouseName).Behav_A50_D10;
        simI_DistH(num_seg)=0;
        simI_Compass(num_seg)=0;
        simI_Behav(num_seg)=0;
        
        for i=1:num_seg
            simI_DistH(i)=dot(NeuP(seg(i),:),DistanceH)/(norm(NeuP(seg(i),:)*norm(DistanceH)));
            simI_Speed(i)=dot(NeuP(seg(i),:),Speed)/(norm(NeuP(seg(i),:)*norm(Speed)));
            simI_Compass(i)=dot(NeuP(seg(i),:),Compass)/(norm(NeuP(seg(i),:)*norm(Compass)));
            simI_Behav(i)=dot(NeuP(seg(i),:),Behav_A50_D10)/(norm(NeuP(seg(i),:)*norm(Behav_A50_D10)));
        end
       
        % variables for shuffled behavioral vectors
        
        shufDistH(orderingSize,frame)=0;
        shufCompass(orderingSize,frame)=0;
        shufBehav(orderingSize,frame)=0;
        
        % creating shuffled behavioral vectors
        for o=1:8
            tempDistH.piece(o).seq=DistanceH((notes(o)+1):notes(o+1));
            tempCompass.piece(o).seq=Compass((notes(o)+1):notes(o+1));
            tempBehav.piece(o).seq=Behav_A50_D10((notes(o)+1):notes(o+1));
        end
        
        for k=1:orderingSize
            shufDistH(k,:)=[[tempDistH.piece(allorders(randorders(k),1)).seq];[tempDistH.piece(allorders(randorders(k),2)).seq];[tempDistH.piece(allorders(randorders(k),3)).seq];[tempDistH.piece(allorders(randorders(k),4)).seq];[tempDistH.piece(allorders(randorders(k),5)).seq];[tempDistH.piece(allorders(randorders(k),6)).seq];[tempDistH.piece(allorders(randorders(k),7)).seq];[tempDistH.piece(allorders(randorders(k),8)).seq]];
            shufCompass(k,:)=[[tempCompass.piece(allorders(randorders(k),1)).seq];[tempCompass.piece(allorders(randorders(k),2)).seq];[tempCompass.piece(allorders(randorders(k),3)).seq];[tempCompass.piece(allorders(randorders(k),4)).seq];[tempCompass.piece(allorders(randorders(k),5)).seq];[tempCompass.piece(allorders(randorders(k),6)).seq];[tempCompass.piece(allorders(randorders(k),7)).seq];[tempCompass.piece(allorders(randorders(k),8)).seq]];
            shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
        end
        
        % similarity index of randomized data: dot sum of two vectors divided by the product of their
        % norm values (5000 randomized behavioral vectors).
        
        simI_shufDistH(num_seg,orderingSize)=0;
        simI_shufCompass(num_seg,orderingSize)=0;
        simI_shufBehav(num_seg,orderingSize)=0;
        
        for i=1:num_seg
            for j=1:orderingSize
                simI_shufDistH(i,j)=dot(NeuP(seg(i),:),shufDistH(j,:))/(norm(NeuP(seg(i),:)*norm(shufDistH(j,:))));
                simI_shufCompass(i,j)=dot(NeuP(seg(i),:),shufCompass(j,:))/(norm(NeuP(seg(i),:)*norm(shufCompass(j,:))));
                simI_shufBehav(i,j)=dot(NeuP(seg(i),:),shufBehav(j,:))/(norm(NeuP(seg(i),:)*norm(shufBehav(j,:))));
            end
        end
        
        % Counting how many randomized data are smaller than actual data
        
        SSeg_DistH(num_seg)=0;
        SSeg_Speed(num_seg)=0;
        SSeg_Compass(num_seg)=0;
        SSeg_Behav(num_seg)=0;
        
        for m=1:num_seg
            SSeg_DistH(m)=numel(find(simI_DistH(m)>simI_shufDistH(m,:)));
            SSeg_Speed(m)=numel(find(simI_Speed(m)>simI_shufSpeed(m,:)));
            SSeg_Compass(m)=numel(find(simI_Compass(m)>simI_shufCompass(m,:)));
            SSeg_Behav(m)=numel(find(simI_Behav(m)>simI_shufBehav(m,:)));
        end

        CSI.("all")=seg(find(SSeg_DistH>4750|SSeg_DistH<250));
        CSI.("allPercent")=numel(CSI.all)/num_seg;
        CSI.("Speed")=seg(find(SSeg_Speed>4750));
        CSI.("Compass")=seg(find(SSeg_Compass>4750));
        CSI.("Behav")=seg(find(SSeg_Behav>4750));
    
        CSI.("allmin_CSI_DistH")=simI_DistH;
        CSI.("allmin_CSI_Behav")=simI_Behav;
        
        % placecell information    
        place.("All_MSI")=MSIseg;
        place.("mAll_MSI")=mean(place.All_MSI(isfinite(place.All_MSI)));
        
        Behav=[];
        for i=1:length(CSI.Behav)
            Behav=[Behav find(seg==CSI.Behav(i))];
        end
        place.("Behav_MSI")=MSIseg(Behav,1);
        place.("mBehav_MSI")=mean(place.Behav_MSI(isfinite(place.Behav_MSI)));
        csiName = "CSI" + suffix;
        placeName = "place" + suffix;
        listCSI.(csiName) = CSI;
        places.(placeName) = place;
    end
    save(filename, "-struct", "listCSI", "-append");
    save(filename, "-struct", "places", "-append");
end