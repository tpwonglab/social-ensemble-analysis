function ensemble_plot(filename)
    if isempty(filename)
            quit(2);
    end
    data = load(filename);
    ensb_Behav = data.ensb_Behav;
    ensb_Behav40 = data.ensb_Behav40;
    NeuP = data.NeuP;

    % Compute mean CaImg data during different close events
    Close = ensb_Behav.CloseEvents;
    NumClose = numel(Close);
    meanNeuP(NumClose, 41) = 0;
    sumNeuP(NumClose) = 0;

    for i=1:NumClose
        meanNeuP(i,:) = sum(NeuP(ensb_Behav40.ensbID,Close(i)-25:Close(i)+15))/numel(ensb_Behav40.ensbID);
        sumNeuP(i) = sum(meanNeuP(i,:));
    end

    mmNeuP = mean(meanNeuP);
    sdNeuP = std(meanNeuP);
    semNeuP = sdNeuP/sqrt(NumClose - 1);

    totalNeuP = sum(mmNeuP);

    % Plot the result
    figure
    plot(mmNeuP)
    hold on
    plot(mmNeuP + semNeuP)
    plot(mmNeuP - semNeuP)
    save(filename, "NumClose", "meanNeuP", "sumNeuP", "mmNeuP", "sdNeuP", "semNeuP", "totalNeuP", "-append")
end