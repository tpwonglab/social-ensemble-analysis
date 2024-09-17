% Input variables in structure array segment
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

% calculate the direction of head vs the enclosure using ispolysw from
% the imaging mapping tools; enclosure on the left of the neck-nose vector
% produces -ve angle; enclosure on the right of the neck-nose vector
% produces +ve angle.
% clc;

Xbox=mean(coordinates.xbox);
Ybox=mean(coordinates.ybox);

[frame,items]=size(coordinates.xn);
% calculate distance
for i=1:frame
coordinates.h2boxDist(i)=sqrt(((coordinates.yh(i)-Ybox)^2)+((coordinates.xh(i)-Xbox)^2));
coordinates.b2boxDist(i)=sqrt(((coordinates.yb(i)-Ybox)^2)+((coordinates.xb(i)-Xbox)^2));
end

% calculate angle vs the box (angle; nose of observer to the box, empty or with a CD1)
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

% calculate compass angle (direction) for the observer mouse
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
else
    coordinates.direction1(i)=d;
end
end


%binning Calmg 2mice downsize
CaImgTime=timestamp.sysClock(timestamp.camNum==segment.CaImgChannel);
x=1;
for i=1:3:numel(CaImgTime)
Catime(x,1)=CaImgTime(i,1);x=x+1;
end
CaImgTime=Catime;

% Input data about characteristics of recording

% prompt = 'Input the # of frames per bin (binwidth): ';
% o = input(prompt);
o = 6;
% Calculate the coordinate of the enclosure
% Bin data according to binwidth
binwidth=(o*(1/30))*1000;
% prompt = 'Input maxtime (last CaImgRawtime): ';
% MaxTime = input(prompt);
MaxTime = segment.CaImgRawtime;
% prompt = 'Input starting frame #: ';
% NeuStart = input(prompt);
NeuStart = segment.NeuStart;
% prompt = 'Input ending frame #: ';
% NeuEnd = input(prompt);
NeuEnd = segment.NeuEnd;
% prompt = 'Input frame number: ';
% CaImgRawFN = input(prompt);
CaImgRawFN = NeuEnd-NeuStart+1;
% prompt = 'Input segment #: ';
% seg = input(prompt);
seg = segment.seg;

binsize=(fix(MaxTime/binwidth)+1);
edgeTime=0:binwidth:(binwidth*binsize);
% CaImgTime=timestamp.sysClock(timestamp.camNum==segment.CaImgChannel);
binnedCaImg=discretize(CaImgTime,edgeTime);

% binning calcium imaging data from coordinates
tempC=zeros(seg,CaImgRawFN);
tempS=zeros(seg,CaImgRawFN);
for j=1:seg
tempC(j,NeuStart:NeuEnd)=neuron.C(j,NeuStart:NeuEnd);
% [pks,locs]=findpeaks(neuron.S(j,NeuStart:NeuEnd));
% tempSpeaks=neuron.S(1,NeuStart:NeuEnd)*0;
% tempSpeaks(1,locs)=1;
tempS(j,NeuStart:NeuEnd)=neuron.S(j,NeuStart:NeuEnd);
end
for i=1:seg
    nc=tempC(i,:).';
    ns=tempS(i,:).';
    tempNeuC(i,:)=accumarray(binnedCaImg,nc,[],@mean);
    tempNeuS(i,:)=accumarray(binnedCaImg,ns,[],@sum);
end
% trimming data to remove data points without CNMF-E analysis

NeuC=tempNeuC(:,binnedCaImg(NeuStart,1):binnedCaImg(NeuEnd,1));
NeuS=tempNeuS(:,binnedCaImg(NeuStart,1):binnedCaImg(NeuEnd,1));

% estimate normNeuC

[seg,frames]=size(NeuC);

normNeuC=[];
for x=1:seg
normNeuC(x,:)=NeuC(x,:)/max(NeuC(x,:));
end

% find active frames > Rmean (Robust mean =
% mean+2SD)

m2sd_frames=fix(frames-0.025*frames);

% Find frames with 

Rmean=[];
activeNeuC=zeros(seg,frames);

Rmean(numel(seg))=0;

for i=1:seg
    tempsortedNeuC=sort(normNeuC(i,:));
    Rmean(i)=mean(tempsortedNeuC(1:m2sd_frames));
    tempactive=find(normNeuC(i,:)>Rmean(i));
    activeNeuC(i,tempactive)=1;
end

% binning behavioral data from coordinates
CaImgTime=timestamp.sysClock(timestamp.camNum==segment.CaImgChannel);
binnedCaImg=discretize(CaImgTime,edgeTime);
BehavTime=timestamp.sysClock(timestamp.camNum==segment.BehavChannel);
binnedBehav=discretize(BehavTime,edgeTime);
tempHeadX=accumarray(binnedBehav,coordinates.xn,[],@mean);
tempHeadY=accumarray(binnedBehav,coordinates.yn,[],@mean);
tempNoseX=accumarray(binnedBehav,coordinates.xh,[],@mean);
tempNoseY=accumarray(binnedBehav,coordinates.yh,[],@mean);
tempAngle=accumarray(binnedBehav,coordinates.angle,[],@mean);
tempDirection=accumarray(binnedBehav,coordinates.direction1,[],@mean);
tempDistanceH=accumarray(binnedBehav,coordinates.h2boxDist,[],@mean);
tempDistanceB=accumarray(binnedBehav,coordinates.b2boxDist,[],@mean);

HeadX=tempHeadX(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
HeadY=tempHeadY(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
NoseX=tempNoseX(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
NoseY=tempNoseY(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
Angle=tempAngle(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
Direction=tempDirection(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
DistanceH=tempDistanceH(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);
DistanceB=tempDistanceB(binnedBehav(segment.BehavStartFN,1):binnedBehav(segment.BehavEndFN,1),1);

View(numel(Direction))=0;

for i=1:numel(Direction)
if Direction(i)>180
View(i)=360-Direction(i);
else
View(i)=Direction(i);
end
end
View=View.';

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

if length(Angle)>length(NeuC)
    Angle(length(Angle))=[];
    HeadX(length(Angle))=[];
    HeadY(length(Angle))=[];
    NoseX(length(Angle))=[];
    NoseY(length(Angle))=[];
    Direction(length(Angle))=[];
    DistanceH(length(Angle))=[];
    DistanceB(length(Angle))=[];
    View(length(Angle))=[];
    Speed(length(Angle))=[];
    Displacement(length(Angle))=[];
else
end

% estimate Close and Far

A50_D10=find(Angle>-50&Angle<50&DistanceH<10);
A50_D25=find(Angle>-50&Angle<50&DistanceH>10&DistanceH<25);
A50_Dfar=find(Angle>-50&Angle<50&DistanceH>25);
A130_D10=find(Angle<-130&DistanceH<10|Angle>130&DistanceH<10);
A130_D25=find(Angle<-130&DistanceH>10&DistanceH<25|Angle>130&DistanceH>10&DistanceH<25);
A130_Dfar=find(Angle<-130&DistanceH>25|Angle>130&DistanceH>25);

frames=length(Angle);

Behav_zeros=zeros(1,frames);
Behav_A50_D10=Behav_zeros;
Behav_A50_D10(A50_D10)=1;
Behav_A50_D25=Behav_zeros;
Behav_A50_D25(A50_D25)=1;
Behav_A50_Dfar=Behav_zeros;
Behav_A50_Dfar(A50_Dfar)=1;
Behav_A130_D10=Behav_zeros;
Behav_A130_D10(A130_D10)=1;
Behav_A130_D25=Behav_zeros;
Behav_A130_D25(A130_D25)=1;
Behav_A130_Dfar=Behav_zeros;
Behav_A130_Dfar(A130_Dfar)=1;

Close=find(DistanceH<10);
Far=find(DistanceH>25);

% save data
segment.binwidth=binwidth;
filename=horzcat(segment.mouseID,'_',segment.session,'.mat');
save(filename,'coordinates','timestamp','segment','HeadX','HeadY','NoseX','NoseY','Xbox','Ybox','Angle','Direction','DistanceH','DistanceB','NeuC','normNeuC','NeuS','Behav_A50_D10','Behav_A50_D25','Behav_A50_Dfar','Behav_A130_D10','Behav_A130_D25','Behav_A130_Dfar','Close','Far','activeNeuC','Rmean','Displacement','Speed','View','neuron');