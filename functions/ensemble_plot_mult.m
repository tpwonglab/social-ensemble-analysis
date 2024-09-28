function ensemble_plot(outputFilename)
    if isempty(outputFilename)
            quit(2);
    end
    data = load(outputFilename);
    segment = data.segment;
    ensb_Behav_A = data.ensb_Behav_A;
    ensb_Behav_N = data.ensb_Behav_N;
    ensb_Behav40_A = data.ensb_Behav40_A;
    ensb_Behav40_N = data.ensb_Behav40_N;
    NeuP = data.NeuP;

    CloseA = ensb_Behav_A.CloseEvents;
    CloseN = ensb_Behav_N.CloseEvents;
    
    NumCloseA=numel(CloseA);
    NumCloseN=numel(CloseN);

    mNeuP_AA(NumCloseA,41)=0;
    mNeuP_NA(NumCloseA,41)=0;
    mNeuP_AN(NumCloseN,41)=0;
    mNeuP_NN(NumCloseN,41)=0;
    sumNeuP_AA(NumCloseA)=0;
    sumNeuP_NA(NumCloseA)=0;
    sumNeuP_AN(NumCloseN)=0;
    sumNeuP_NN(NumCloseN)=0;
    
    for i=1:NumCloseA
        mNeuP_AA(i,:)=sum(NeuP(ensb_Behav40_A.ensbID,CloseA(i)-25:CloseA(i)+15))/numel(ensb_Behav40_A.ensbID);
        mNeuP_NA(i,:)=sum(NeuP(ensb_Behav40_N.ensbID,CloseA(i)-25:CloseA(i)+15))/numel(ensb_Behav40_N.ensbID);
        sumNeuP_AA(i)=sum(mNeuP_AA(i,:));
        sumNeuP_NA(i)=sum(mNeuP_NA(i,:));
    end
    
    mmNeuP_AA=mean(mNeuP_AA);
    mmNeuP_NA=mean(mNeuP_NA);
    sdNeuP_AA=std(mNeuP_AA);
    sdNeuP_NA=std(mNeuP_NA);
    semNeuP_AA=sdNeuP_AA/sqrt(NumCloseA-1);
    semNeuP_NA=sdNeuP_NA/sqrt(NumCloseA-1);
    
    for j=1:NumCloseN
        mNeuP_AN(j,:)=sum(NeuP(ensb_Behav40_A.ensbID,CloseN(j)-25:CloseN(j)+15))/numel(ensb_Behav40_A.ensbID);
        mNeuP_NN(j,:)=sum(NeuP(ensb_Behav40_N.ensbID,CloseN(j)-25:CloseN(j)+15))/numel(ensb_Behav40_N.ensbID);
        sumNeuP_AN(j)=sum(mNeuP_AN(j,:));
        sumNeuP_NN(j)=sum(mNeuP_NN(j,:));
    end
    
    mmNeuP_AN=mean(mNeuP_AN);
    mmNeuP_NN=mean(mNeuP_NN);
    sdNeuP_AN=std(mNeuP_AN);
    sdNeuP_NN=std(mNeuP_NN);
    semNeuP_AN=sdNeuP_AN/sqrt(NumCloseN-1);
    semNeuP_NN=sdNeuP_NN/sqrt(NumCloseN-1);
    
    totalNeuP_AA=sum(mmNeuP_AA);
    totalNeuP_NA=sum(mmNeuP_NA);
    totalNeuP_NN=sum(mmNeuP_NN);
    totalNeuP_AN=sum(mmNeuP_AN);
    
    [hA,pA]=ttest(sumNeuP_AA,sumNeuP_NA);
    [hN,pN]=ttest(sumNeuP_NN,sumNeuP_AN);

    figure
    plot(mmNeuP_AA)
    hold on
    plot(mmNeuP_AA+semNeuP_AA)
    plot(mmNeuP_AA-semNeuP_AA)
    plot(mmNeuP_AN)
    plot(mmNeuP_AN+semNeuP_AN)
    plot(mmNeuP_AN-semNeuP_AN)
    % filename=horzcat(segment.mouseID,'_',segment.session,'_AAvAN40.fig');
    % saveas(gcf,filename)
    
    figure
    plot(mmNeuP_NN)
    hold on
    plot(mmNeuP_NN+semNeuP_NN)
    plot(mmNeuP_NN-semNeuP_NN)
    plot(mmNeuP_NA)
    plot(mmNeuP_NA+semNeuP_NA)
    plot(mmNeuP_NA-semNeuP_NA)
    % filename=horzcat(segment.mouseID,'_',segment.session,'_freqAnalysis40.mat');

    save(outputFilename,"NumCloseA","NumCloseN","NumCloseA","NumCloseN", "-append");
    save(outputFilename,"mNeuP_AA","mNeuP_AN","mNeuP_NA","mNeuP_NN","-append");
    save(outputFilename,"mmNeuP_AA","mmNeuP_AN","mmNeuP_NA","mmNeuP_NN","-append");
    save(outputFilename,"sdNeuP_AA","sdNeuP_AN","sdNeuP_NA","sdNeuP_NN","-append");
    save(outputFilename,"semNeuP_AA","semNeuP_AN","semNeuP_NA","semNeuP_NN","-append");
    save(outputFilename,"totalNeuP_AA","totalNeuP_AN","totalNeuP_NA","totalNeuP_NN","-append");
    save(outputFilename,"hA","hN","pA","pN","-append");
end