function outputFilename = binning(segmentFilename, dataFilename, isSingleMouse)
    %% Define segment
    if isempty(segmentFilename)
        disp("Create new segment file")
        segment = askSegmentPrompt;
    else
        disp("Create new segment file")
        disp("Load " + segmentFilename + " file")
        segment = load(segmentFilename);
    end
    if isempty(dataFilename)
        disp("Missing data filename");
        quit(2);
    end

    %% Load initial raw data Calcium Imaging and Real-time Video
    data = load(dataFilename);
    coordinates = data.coordinates;
    neuron = data.neuron;
    timestamp = data.timestamp;

    %% Bin Calcium Imaging Video
    % Sort timestamps only for Calcium imaging timestamps
    CaImgTime=timestamp.sysClock(timestamp.camNum==segment.CaImgChannel);
    x=1;
    for i=1:3:numel(CaImgTime)
        Catime(x,1)=CaImgTime(i,1);x=x+1;
    end
    CaImgTime=Catime;

    % Calculate the coordinate of the enclosure; bin data according to binwidth
    o = 6;
    binwidth=(o*(1/30))*1000;
    MaxTime = segment.CaImgRawtime;
    NeuStart = segment.NeuStart;
    NeuEnd = segment.NeuEnd;
    CaImgRawFN = NeuEnd-NeuStart+1;
    seg = segment.seg;
    binsize=(fix(MaxTime/binwidth)+1);
    edgeTime=0:binwidth:(binwidth*binsize);
    binnedCaImg=discretize(CaImgTime,edgeTime);
    % Binning calcium imaging data from coordinates
    tempS=zeros(seg,CaImgRawFN);
    for j=1:seg
        tempS(j,NeuStart:NeuEnd)=neuron.S(j,NeuStart:NeuEnd);
    end
    for i=1:seg
        ns=tempS(i,:).';
        tempNeuS(i,:)=accumarray(binnedCaImg,ns,[],@sum);
    end
    % trimming data to remove data points without CNMF-E analysis
    NeuS=tempNeuS(:,binnedCaImg(1,1):binnedCaImg(CaImgRawFN,1));

    %% Bin Behavioral Recording
    if isSingleMouse == 1
        [frame, items]=size(coordinates.xn);
        Xbox=mean(coordinates.xbox);
        Ybox=mean(coordinates.ybox);
        % Calculate the head distance between mouse and box
        for i=1:frame
            coordinates.h2boxDist(i)=sqrt(((coordinates.yh(i)-Ybox)^2)+((coordinates.xh(i)-Xbox)^2));
        end
        % Calculate the head angle between mouse and box
        for i=1:frame
            serX=[Xbox coordinates.xh(i) coordinates.xn(i)];
            serY=[Ybox coordinates.yh(i) coordinates.yn(i)];
            c=sqrt((Ybox-coordinates.yn(i))^2+(Xbox-coordinates.xn(i))^2);
            a=sqrt((Ybox-coordinates.yh(i))^2+(Xbox-coordinates.xh(i))^2);
            b=sqrt((coordinates.yh(i)-coordinates.yn(i))^2+(coordinates.xh(i)-coordinates.xn(i))^2);
            d=rad2deg(acos((b^2+c^2-a^2)/(2*b*c)));
            tf=ispolycw(serX,serY);
            if tf>0
                coordinates.angle(i)=d;
            else
                coordinates.angle(i)=0-d;
            end
        end
        % Calculate the compass angle (direction) for the observer mouse
        Compass(numel(coordinates.xn)) = 0; Compass = Compass.';
        for i=1:frame
            serX=[coordinates.xn(i) coordinates.xh(i) coordinates.xn(i)];
            serY=[0 coordinates.yh(i) coordinates.yn(i)];
            c=sqrt((0-coordinates.yn(i))^2+(coordinates.xn(i)-coordinates.xn(i))^2);
            a=sqrt((0-coordinates.yh(i))^2+(coordinates.xn(i)-coordinates.xh(i))^2);
            b=sqrt((coordinates.yh(i)-coordinates.yn(i))^2+(coordinates.xh(i)-coordinates.xn(i))^2);
            d=rad2deg(acos((b^2+c^2-a^2)/(2*b*c)));
            tf=ispolycw(serX,serY);
            if tf>0
                coordinates.direction1(i)=360-d;
                Compass(i) = 360 - d;
            else
                coordinates.direction1(i)=d;
                Compass(i) = d;
            end
        end

        BehavTime=timestamp.sysClock(timestamp.camNum==segment.BehavChannel);
        binnedBehav=discretize(BehavTime,edgeTime);
        tempHeadX=accumarray(binnedBehav,coordinates.xn,[],@mean);
        tempHeadY=accumarray(binnedBehav,coordinates.yn,[],@mean);
        tempNoseX=accumarray(binnedBehav,coordinates.xh,[],@mean);
        tempNoseY=accumarray(binnedBehav,coordinates.yh,[],@mean);
        tempAngle=accumarray(binnedBehav,coordinates.angle,[],@mean);
        tempDistanceH=accumarray(binnedBehav,coordinates.h2boxDist,[],@mean);

        HeadX=tempHeadX(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
        HeadY=tempHeadY(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
        NoseX=tempNoseX(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
        NoseY=tempNoseY(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
        Angle=tempAngle(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
        DistanceH=tempDistanceH(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);

        % Speed and Displacement
        frames=numel(Angle);
        Speed(frames)=0;
        Displacement(frames)=0;
        for i=1:frames-1
            Speed(i+1)=sqrt((HeadX(i+1)-HeadX(i))^2+(HeadY(i+1)-HeadY(i))^2)*5;
            Displacement(i+1)=sqrt((HeadX(i+1)-HeadX(i))^2+(HeadY(i+1)-HeadY(i))^2);
        end
        Speed=Speed.';
        Displacement=Displacement.';

        % Remove extra frames
        if length(Angle)>length(NeuS)
            Angle(length(Angle))=[];
            HeadX(length(Angle))=[];
            HeadY(length(Angle))=[];
            NoseX(length(Angle))=[];
            NoseY(length(Angle))=[];
            DistanceH(length(Angle))=[];
            Speed(length(Angle))=[];
            Displacement(length(Angle))=[];
        end

        % Generate NeuP from findpeaks
        NeuP(seg,length(HeadX))=0;
        for i=1:seg
            [peaks,locs]=findpeaks(NeuS(i,:));
            for j=1:numel(peaks)
                NeuP(i,locs(j))=1;
            end
        end

        % estimate Close and Far
        A50_D10=find(Angle>-50&Angle<50&DistanceH<10);
        A50_Dfar=find(Angle>-50&Angle<50&DistanceH>25);
        frames=length(Angle);
        
        Behav_zeros=zeros(1,frames);
        Behav_A50_D10=Behav_zeros;
        Behav_A50_D10(A50_D10)=1;
        Behav_A50_Dfar=Behav_zeros;
        Behav_A50_Dfar(A50_Dfar)=1;
        segment.binwidth=binwidth;
        filename=segment.mouseID + '_' + segment.session + '.mat';
        if exist(filename, "file")
            delete(filename);
        end
        outputFilename = filename;
        save(filename, 'coordinates', 'timestamp', 'neuron', 'segment', 'HeadX', 'HeadY', 'NoseX', 'NoseY', 'Xbox', 'Ybox', 'Angle', 'Compass', 'DistanceH', 'NeuS', 'NeuP', 'Behav_A50_D10', 'Behav_A50_Dfar', 'Displacement', 'Speed');
        return
    end
    % Binning behavioral data from coordinates
    BehavTime=timestamp.sysClock(timestamp.camNum==segment.BehavChannel);
    binnedBehav=discretize(BehavTime,edgeTime);
    tempHeadX=accumarray(binnedBehav,coordinates.xn,[],@mean);
    tempHeadX2=accumarray(binnedBehav,coordinates.xn2,[],@mean);
    tempHeadY=accumarray(binnedBehav,coordinates.yn,[],@mean);
    tempHeadY2=accumarray(binnedBehav,coordinates.yn2,[],@mean);
    tempNoseX=accumarray(binnedBehav,coordinates.xh,[],@mean);
    tempNoseX2=accumarray(binnedBehav,coordinates.xh2,[],@mean);
    tempNoseY=accumarray(binnedBehav,coordinates.yh,[],@mean);
    tempNoseY2=accumarray(binnedBehav,coordinates.yh2,[],@mean);
    
    HeadX=tempHeadX(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    HeadX2=tempHeadX2(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    HeadY=tempHeadY(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    HeadY2=tempHeadY2(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    NoseX=tempNoseX(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    NoseX2=tempNoseX2(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    NoseY=tempNoseY(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
    NoseY2=tempNoseY2(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);

    % Speed and Displacement
    frames=numel(HeadX);
    Speed(frames)=0;
    Displacement(frames)=0;
    for i=1:frames-1
        Speed(i+1)=sqrt((HeadX(i+1)-HeadX(i))^2+(HeadY(i+1)-HeadY(i))^2)*5;
        Displacement(i+1)=sqrt((HeadX(i+1)-HeadX(i))^2+(HeadY(i+1)-HeadY(i))^2);
    end
    Speed=Speed.';
    Displacement=Displacement.';
    
    % Remove extra frames
    
    if length(HeadX)>length(NeuS)
        HeadX(length(HeadX))=[];
        HeadX2(length(HeadX))=[];
        HeadY(length(HeadX))=[];
        HeadY2(length(HeadX))=[];
        NoseX(length(HeadX))=[];
        NoseX2(length(HeadX))=[];
        NoseY(length(HeadX))=[];
        NoseY2(length(HeadX))=[];
        Speed(length(HeadX))=[];
        Displacement(length(HeadX))=[];
    end
    
    % Generate NeuP from findpeaks
    NeuP(seg,length(HeadX))=0;
    for i=1:seg
        [peaks,locs]=findpeaks(NeuS(i,:));
        for j=1:numel(peaks)
            NeuP(i,locs(j))=1;
        end
    end
    
    % calculate DistanceH, Angle (angle to the head of mouse 2), Compass (compass angle, 0 at
    % North, +ve clockwise, -ve anticlockwise)
    
    DistanceH(numel(HeadX))=0;DistanceH=DistanceH.';
    Angle(numel(HeadX))=0;Angle=Angle.';Angle2=Angle;
    Compass(numel(HeadX))=0;Compass=Compass.';
    for i=1:numel(HeadX)
        DistanceH(i)=sqrt((HeadX(i)-HeadX2(i))^2+(HeadY(i)-HeadY2(i))^2);
    
        a=sqrt((NoseX(i)-HeadX2(i))^2+(NoseY(i)-HeadY2(i))^2);
        b=sqrt((NoseX(i)-HeadX(i))^2+(NoseY(i)-HeadY(i))^2);
        c=sqrt((HeadX(i)-HeadX2(i))^2+(HeadY(i)-HeadY2(i))^2);
        Angle(i)=rad2deg(acos((b^2+c^2-a^2)/(2*b*c)));
    
        a=sqrt((NoseX2(i)-HeadX(i))^2+(NoseY2(i)-HeadY(i))^2);
        b=sqrt((NoseX2(i)-HeadX2(i))^2+(NoseY2(i)-HeadY2(i))^2);
        c=sqrt((HeadX2(i)-HeadX(i))^2+(HeadY2(i)-HeadY(i))^2);
        Angle2(i)=rad2deg(acos((b^2+c^2-a^2)/(2*b*c)));
    
        serX=[NoseX(i) HeadX(i) HeadX(i)];
        serY=[NoseY(i) HeadY(i) (NoseY(i)+HeadY(i))];
        e=sqrt((NoseX(i)-HeadX(i))^2+(NoseY(i)-(NoseY(i)+HeadY(i)))^2);
        f=sqrt((NoseX(i)-HeadX(i))^2+(NoseY(i)-HeadY(i))^2);
        g=sqrt((HeadX(i)-HeadX(i))^2+(HeadY(i)-(NoseY(i)+HeadY(i)))^2);
        h=rad2deg(acos((f^2+g^2-e^2)/(2*f*g)));
        tf=ispolycw(serX,serY);
        if tf>0
            Compass(i)=h;
        else
            Compass(i)=0-h;
        end
    end
    A50_D10=find(Angle<50&DistanceH<10);
    A50_Dfar=find(Angle<50&DistanceH>25);
    
    frames=length(Angle);
    
    Behav_zeros=zeros(1,frames);
    Behav_A50_D10=Behav_zeros;
    Behav_A50_D10(A50_D10)=1;Behav_A50_D10=Behav_A50_D10.';
    
    Behav_A50_Dfar=Behav_zeros;
    Behav_A50_Dfar(A50_Dfar)=1;Behav_A50_Dfar=Behav_A50_Dfar.';
    segment.binwidth=binwidth;
    filename=segment.mouseID + '_' + segment.session + '.mat';
    if exist(filename, "file")
        delete(filename);
    end
    outputFilename = filename;
    save(filename,'coordinates','timestamp', 'neuron','segment','HeadX','HeadX2','HeadY','HeadY2','NoseX','NoseX2','NoseY','NoseY2','Angle','Angle2','Compass','DistanceH','NeuS','NeuP','Behav_A50_D10','Behav_A50_Dfar','Displacement','Speed');
end

function segment = askSegmentPrompt()
    prompt='Please input the Mouse ID:';
    segment.mouseID = input(prompt,'s');
    prompt='Please input the session name (Def2, SIA, OF, etc):';
    segment.session = input(prompt,'s');
    prompt='Please input the camera channel for CaImg:';
    segment.CaImgChannel = input(prompt);
    prompt='Please input the camera channel for Behav:';
    segment.BehavChannel = input(prompt);
    prompt='Please input the raw # of frames of CaImg before CNMF-E:';
    segment.CaImgRawFN = input(prompt);
    prompt='Please input the total time (in ms) of all frames of CaImg before CNMF-E:';
    segment.CaImgRawtime = input(prompt);
    prompt='Please input the raw # of frames of Behav:';
    segment.BehavRawFN = input(prompt);
    prompt='Please input the total time (in ms) of all frames of Behav:';
    segment.BehavRawtime = input(prompt);
    prompt='Please input the starting frame # of CaImg that was analyzed by CNMF-E:';
    segment.CaImgStartFN = input(prompt);
    prompt='Please input the starting frame time (in ms) of CaImg that was analyzed by CNMF-E:';
    segment.CaImgStartFtime = input(prompt);
    prompt='Please input the ending frame # of CaImg that was analyzed by CNMF-E:';
    segment.CaImgEndFN = input(prompt);
    prompt='Please input the starting frame # of Behav with analyzed CaImg data:';
    segment.BehavStartFN = input(prompt);
    prompt='Please input the starting frame time (in ms) of Behav with analyzed CaImg data:';
    segment.BehavStartFtime = input(prompt);
    prompt='Please input the ending frame # of Behav with analyzed CaImg data:';
    segment.BehavEndFN = input(prompt);
    prompt='Please input the segments:';
    segment.seg = input(prompt);
    prompt='Please input the starting frame of Neuron.S:';
    segment.NeuStart = input(prompt);
    prompt='Please input the ending frame of Neuron.S:';
    segment.NeuEnd = input(prompt);
end