function ensemble(segmentFilename, dataFilename, experimentName)
    if isempty(dataFilename)
            quit(2);
    end
    data = load(dataFilename);
    shuffleRuns = 10000;
    suffixes = [""];
    if experimentName == "SDT"
        suffixes = ["_A", "_N"];
    end
    for k=1:length(suffixes)
        clear ensb_Behav ensb_Behav10 ensb_Behav20 ensb_Behav30 ensb_Behav40 ensb_Behav50
        suffix = suffixes(k);
        mouseName = "mouse" + suffix;
        Behav_A50_D10 = data.(mouseName).Behav_A50_D10;
        Behav_A50_Dfar = data.(mouseName).Behav_A50_Dfar;
        if experimentName == "SI1"
            load(segmentFilename);
            Behav_A50_D10(startS2:endS2)=0;
            Behav_A50_Dfar(startS2:endS2)=0;
        elseif experimentName == "SI2"
            load(segmentFilename);
            Behav_A50_D10(1:startS2-1)=0;
            Behav_A50_Dfar(1:startS2-1)=0;
        end
        NeuP = data.NeuP;
        DistanceH = data.(mouseName).DistanceH;
        Compass = data.(mouseName).Compass;
        Speed = data.Speed;
        Angle = data.(mouseName).Angle;
    
        counts=Behav_A50_D10;
        counts(1:25,1)=0;
        counts((length(counts)-15:length(counts)))=0;
        counts(1:1,1)=0;

        counts(length(counts):length(counts), 1) = 0;
        
        Behav_A50_D10=counts;
        
        for i=1:(length(Behav_A50_D10)-1)
            if counts(i)>0&&counts(i+1)>0
                counts(i+1)=counts(i)+counts(i+1);
            end
        end
        
        tempAns=input("Initiate ensemble analysis (Yes=1; No=0)? ");

        if tempAns==1

        CloseEvents=find(counts==5);
        CloseEventsFrames(length(Behav_A50_D10))=0;
        for i=1:numel(CloseEvents)
            CloseEventsFrames((CloseEvents(i)-10):CloseEvents(i))=1;
        end
        
        [seg,frame]=size(NeuP);
        seg = 1:seg;
        num_seg = numel(seg);
        activeClose(num_seg,numel(CloseEvents))=0;
        
        for i=1:num_seg
            for j=1:numel(CloseEvents)
                tempsum=sum(NeuP(i,(CloseEvents(j)-25):(CloseEvents(j)+15)));
                if tempsum>0
                    activeClose(i,j)=1;
                end
            end
        end
        
        activeClose_seg=sum(activeClose,2);
        activeClose_SI=sum(activeClose,1);
        
        ct=1;
        
        x=[];
        
        DistanceH=50-DistanceH;
        
        Iter=[];

        % Can only run for 5 times for now (but if possible, run 9 times)
        Iter = input("Enter number of iterations: ");
        for x=0.1:0.1:Iter/10
        
            result.percent(ct)=x;
            resultSpec.percent=x;
        
            ensbsize=fix(numel(activeClose_SI)*x);
            ensbTemp=find(activeClose_seg>ensbsize);
            ensb_result=intersect(ensbTemp,seg);
        
            ensb_NeuP=sum(NeuP(ensb_result,:))/numel(ensb_result);
            simI_ensbAngle=dot(ensb_NeuP,Angle)/(norm(ensb_NeuP)*norm(Angle));
            simI_ensbBehav_A50_D10=dot(ensb_NeuP,Behav_A50_D10)/(norm(ensb_NeuP)*norm(Behav_A50_D10));
            simI_ensbBehav_A50_Dfar=dot(ensb_NeuP,Behav_A50_Dfar)/(norm(ensb_NeuP)*norm(Behav_A50_Dfar));
        
            simI_ensbCompass=dot(ensb_NeuP,Compass)/(norm(ensb_NeuP)*norm(Compass));
            simI_ensbSpeed=dot(ensb_NeuP,Speed)/(norm(ensb_NeuP)*norm(Speed));
            simI_ensbDistanceH=dot(ensb_NeuP,DistanceH)/(norm(ensb_NeuP)*norm(DistanceH));
        
        
            % Angle
            
            % shuffle behavioral vectors into 8 pieces
            order=1:8;
            allorders=perms(order);
            [items,ss]=size(allorders);
            randorders=randperm(items,shuffleRuns);
            note_interval=fix(frame/8);
            notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Behav_A50_D10)];
            
            % variables for shuffled behavioral vectors
            
            shufBehav(shuffleRuns,frame)=0;
            
            % creating shuffled behavioral vectors
            
            for o=1:8
                tempBehav.piece(o).seq=Angle((notes(o)+1):notes(o+1));
            end
            
            for k=1:shuffleRuns
                shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
            end
            
            % similarity index of randomized data: dot sum of two vectors divided by the product of their
            % norm values (10000 randomized behavioral vectors).
            
            simI_shufBehav_ensb(shuffleRuns)=0;
            
            for i=1:shuffleRuns
                simI_shufBehav_ensb(i)=dot(ensb_NeuP,shufBehav(i,:))/(norm(ensb_NeuP)*norm(shufBehav(i,:)));
            end
            
            result.over_Angle(ct)=numel(find(simI_shufBehav_ensb>simI_ensbAngle));
            resultSpec.over_Angle=result.over_Angle(ct);
        
            % Behav_A50_D10
            
            % shuffle behavioral vectors into 8 pieces
            
            order=1:8;
            allorders=perms(order);
            [items,ss]=size(allorders);
            randorders=randperm(items,shuffleRuns);
            note_interval=fix(frame/8);
            notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Behav_A50_D10)];
            
            % variables for shuffled behavioral vectors
            
            shufBehav(shuffleRuns,frame)=0;
            
            % creating shuffled behavioral vectors
            
            for o=1:8
                tempBehav.piece(o).seq=Behav_A50_D10((notes(o)+1):notes(o+1));
            end
            
            for k=1:shuffleRuns
                shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
            end
            
            % similarity index of randomized data: dot sum of two vectors divided by the product of their
            % norm values (10000 randomized behavioral vectors).
            
            simI_shufBehav_ensb(shuffleRuns)=0;
            
            for i=1:shuffleRuns
                simI_shufBehav_ensb(i)=dot(ensb_NeuP,shufBehav(i,:))/(norm(ensb_NeuP)*norm(shufBehav(i,:)));
            end
            
            result.over_Behav_A50_D10(ct)=numel(find(simI_shufBehav_ensb>simI_ensbBehav_A50_D10));
            resultSpec.over_Behav_A50_D10=result.over_Behav_A50_D10(ct);
                
        
            % Behav_A50_Dfar
            
            % shuffle behavioral vectors into 8 pieces
            
            order=1:8;
            allorders=perms(order);
            [items,ss]=size(allorders);
            randorders=randperm(items,shuffleRuns);
            note_interval=fix(frame/8);
            notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Behav_A50_Dfar)];
            
            % variables for shuffled behavioral vectors
            
            shufBehav(shuffleRuns,frame)=0;
            
            % creating shuffled behavioral vectors
            
            for o=1:8
                tempBehav.piece(o).seq=Behav_A50_Dfar((notes(o)+1):notes(o+1));
            end
            
            for k=1:shuffleRuns;
                shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
            end
            
            % similarity index of randomized data: dot sum of two vectors divided by the product of their
            % norm values (10000 randomized behavioral vectors).
            
            simI_shufBehav_ensb(shuffleRuns)=0;
            
            for i=1:shuffleRuns
                simI_shufBehav_ensb(i)=dot(ensb_NeuP,shufBehav(i,:))/(norm(ensb_NeuP)*norm(shufBehav(i,:)));
            end
            
            result.over_Behav_A50_Dfar(ct)=numel(find(simI_shufBehav_ensb>simI_ensbBehav_A50_Dfar));
            resultSpec.over_Behav_A50_Dfar=result.over_Behav_A50_Dfar(ct);
        
            % Compass
            
            % shuffle behavioral vectors into 8 pieces
            
            order=1:8;
            allorders=perms(order);
            [items,ss]=size(allorders);
            randorders=randperm(items,shuffleRuns);
            note_interval=fix(frame/8);
            notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Compass)];
            
            % variables for shuffled behavioral vectors
            
            shufBehav(shuffleRuns,frame)=0;
            
            % creating shuffled behavioral vectors
            
            for o=1:8
                tempBehav.piece(o).seq=Compass((notes(o)+1):notes(o+1));
            end
            
            for k=1:shuffleRuns;
                shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
            end
            
            % similarity index of randomized data: dot sum of two vectors divided by the product of their
            % norm values (shuffleRuns randomized behavioral vectors).
            
            simI_shufBehav_ensb(shuffleRuns)=0;
            
            for i=1:shuffleRuns
                simI_shufBehav_ensb(i)=dot(ensb_NeuP,shufBehav(i,:))/(norm(ensb_NeuP)*norm(shufBehav(i,:)));
            end
            
            result.over_Compass(ct)=numel(find(simI_shufBehav_ensb>simI_ensbCompass));
            resultSpec.over_Compass=result.over_Compass(ct);
        
            % Speed
            
            % shuffle behavioral vectors into 8 pieces
            % 
            order=1:8;
            allorders=perms(order);
            [items,ss]=size(allorders);
            randorders=randperm(items,shuffleRuns);
            note_interval=fix(frame/8);
            notes=[0 note_interval note_interval*2 note_interval*3 note_interval*4 note_interval*5 note_interval*6 note_interval*7 length(Speed)];
            
            % variables for shuffled behavioral vectors
            
            shufBehav(shuffleRuns,frame)=0;
            
            % creating shuffled behavioral vectors
            
            for o=1:8
                tempBehav.piece(o).seq=Speed((notes(o)+1):notes(o+1));
            end
            
            for k=1:shuffleRuns;
                shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
            end
            
            % similarity index of randomized data: dot sum of two vectors divided by the product of their
            % norm values (10000 randomized behavioral vectors).
            
            simI_shufBehav_ensb(shuffleRuns)=0;
            
            for i=1:shuffleRuns
                simI_shufBehav_ensb(i)=dot(ensb_NeuP,shufBehav(i,:))/(norm(ensb_NeuP)*norm(shufBehav(i,:)));
            end
            
            result.over_Speed(ct)=numel(find(simI_shufBehav_ensb>simI_ensbSpeed));
            resultSpec.over_Speed=result.over_Speed(ct);
        
            % DistanceH
            
            % creating shuffled behavioral vectors: DistanceH
            
            for o=1:8
                tempBehav.piece(o).seq=DistanceH((notes(o)+1):notes(o+1));
            end
            
            for k=1:shuffleRuns;
                shufBehav(k,:)=[[tempBehav.piece(allorders(randorders(k),1)).seq];[tempBehav.piece(allorders(randorders(k),2)).seq];[tempBehav.piece(allorders(randorders(k),3)).seq];[tempBehav.piece(allorders(randorders(k),4)).seq];[tempBehav.piece(allorders(randorders(k),5)).seq];[tempBehav.piece(allorders(randorders(k),6)).seq];[tempBehav.piece(allorders(randorders(k),7)).seq];[tempBehav.piece(allorders(randorders(k),8)).seq]];
            end
            
            % similarity index of randomized data: dot sum of two vectors divided by the product of their
            % norm values (10000 randomized behavioral vectors).
            
            simI_shufBehav_ensb(shuffleRuns)=0;
            
            for i=1:shuffleRuns
                simI_shufBehav_ensb(i)=dot(ensb_NeuP,shufBehav(i,:))/(norm(ensb_NeuP)*norm(shufBehav(i,:)));
            end
            
            result.over_DistanceH(ct)=numel(find(simI_shufBehav_ensb>simI_ensbDistanceH));
            resultSpec.over_DistanceH=result.over_DistanceH(ct);
                
        
            % Organizing results
        
            result.ensb(ct)=numel(ensb_result);
            result.ensbPercent(ct)=numel(ensb_result)/num_seg;
            result.CloseEvents=CloseEvents;
            result.ensbCSI_Angle(ct)=simI_ensbAngle;
            result.ensbCSI_Behav_A50_D10(ct)=simI_ensbBehav_A50_D10;
            result.ensbCSI_Behav_A50_Dfar(ct)=simI_ensbBehav_A50_Dfar;
            result.ensbCSI_Compass(ct)=simI_ensbCompass;
            result.ensbCSI_Speed(ct)=simI_ensbSpeed;
            result.ensbCSI_DistanceH(ct)=simI_ensbDistanceH;
            result.events=i;
            result.sig_events=i/20;
            
            % Calculate for specific ct
            resultSpec.ensb=result.ensb(ct);
            resultSpec.ensbID=ensb_result;
            resultSpec.ensbPercent=result.ensbPercent(ct);
            resultSpec.CloseEvents=CloseEvents;
            resultSpec.ensbCSI_Angle=result.ensbCSI_Angle(ct);
            resultSpec.ensbCSI_Behav_A50_D10=result.ensbCSI_Behav_A50_D10(ct);
            resultSpec.ensbCSI_Behav_A50_Dfar=result.ensbCSI_Behav_A50_Dfar(ct);
            resultSpec.ensbCSI_Compass=result.ensbCSI_Compass(ct);
            resultSpec.ensbCSI_Speed=result.ensbCSI_Speed(ct);
            resultSpec.ensbCSI_DistanceH=result.ensbCSI_DistanceH(ct);
            result.events=i;
            result.sig_events=i/20;
        
            % Save specific ct
            % It may crash, save what you can (suggestion: use try-catch)
            % Save each session at a time (Behav10, Behav20, ...)
            if ct == 1
                ensb_Behav10=resultSpec;
            elseif ct==2
                ensb_Behav20=resultSpec;
            elseif ct==3
                ensb_Behav30=resultSpec;
            elseif ct==4
                ensb_Behav40=resultSpec;
            elseif ct==5
                ensb_Behav50=resultSpec;
            end
        
            ct=ct+1;
            ensb_Behav=result;
        end
        ensb = struct();
        if exist('ensb_Behav10')
            ensb.("ensb_Behav10" + suffix) = struct();
            ensb.("ensb_Behav10" + suffix) = ensb_Behav10;
        end
        if exist('ensb_Behav20')
            ensb.("ensb_Behav20" + suffix) = struct();
            ensb.("ensb_Behav20" + suffix) = ensb_Behav20;
        end
        if exist('ensb_Behav30')
            ensb.("ensb_Behav30" + suffix) = struct();
            ensb.("ensb_Behav30" + suffix) = ensb_Behav30;
        end
        if exist('ensb_Behav40')
            ensb.("ensb_Behav40" + suffix) = struct();
            ensb.("ensb_Behav40" + suffix) = ensb_Behav40;
        end
        if exist('ensb_Behav50')
            ensb.("ensb_Behav50" + suffix) = struct();
            ensb.("ensb_Behav50" + suffix) = ensb_Behav50;
        end
        ensb.("ensb_Behav" + suffix) = struct();
        ensb.("ensb_Behav" + suffix) = ensb_Behav;
        save(dataFilename, "-struct", "ensb", "-append");
        if isempty(extract(experimentName,'SI'))==0&exist('endS2new')
            save(dataFilename, "startS1", "startS2", "endS1", "endS2", "endS2new", "-append")
        elseif isempty(extract(experimentName,'SI'))==0
            save(dataFilename, "startS1", "startS2", "endS1", "endS2", "-append")
        end
    end
end

