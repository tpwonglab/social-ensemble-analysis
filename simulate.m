clear main;
clc;
main()

function main
    disp("Let's get started...");

    % User Inputs
    testingFilename = "CSDS/1144_SDT.mat";
    experimentName = "SDT";
    segmentFilename = "data/init_" + experimentName + "_seg.mat";
    dataFilename = "data/init_" + experimentName + "_data.mat";
    startFrame = 1;

    addpath("functions");

    disp("0. Setup Source2D visualization");
    tic
    setupSource()
    toc

    disp("1. Create initial data");
    tic
    createTestData(testingFilename, segmentFilename, dataFilename, experimentName);
    toc

    disp("2. Create binning between Calcium Imaging and Experiment Video");
    tic
    outputFilename = binning(segmentFilename, dataFilename, experimentName);
    toc

    disp("3. Apply CSI");
    tic
    csi(outputFilename, experimentName);
    toc

    disp("4. Ensemble behavioural");
    tic
    ensemble(outputFilename, startFrame, experimentName);
    toc

    disp("5. Flat map all scenario specific data");
    data = load(outputFilename);
    suffixes = [""];
    if experimentName == "SDT"
        suffixes = ["_A", "_N"];
    end
    tic
    for k = 1:length(suffixes)
        flatten_mouse = struct();
        suffix = suffixes(k);
        mouseName = "mouse" + suffix;
        mouse = data.(mouseName);
        mouseFields = fieldnames(mouse);
        for i = 1:numel(mouseFields)
            name = mouseFields{i};
            flatten_mouse.(name + suffix) = mouse.(name);
        end
        save(outputFilename, "-struct", "flatten_mouse", "-append");
    end
    toc

    disp("6. Remove redundant data from mat file");
    tic
    data = load(outputFilename);
    for k = 1:length(suffixes)
        suffix = suffixes(k);
        fieldName = "mouse" + suffix;
        if isfield(data, fieldName) == 1
            data = rmfield(data, fieldName);
        end
    end
    save(outputFilename, "-struct", "data");
    toc

    disp("7. Plot behavioural data");
    tic
    if experimentName == "SDT"
        ensemble_plot_mult(outputFilename);
    else
        ensemble_plot(outputFilename);
    end
    toc

    disp("Completed.");
end

function setupSource
    run CNMF_E/cnmfe_setup.m;
end

function createTestData(initialData, segmentFilename, dataFilename, experimentName)
    if exist(segmentFilename, "file") && exist(dataFilename, "file")
        return
    end
    data = load(initialData);
    defineInitExp(data.segment, segmentFilename, experimentName);
    defineInitData(data, dataFilename);
end

function defineInitExp(data, filename, experimentName)
    mouseID = "1";
    session = experimentName;
    CaImgChannel = data.CaImgChannel;
    BehavChannel = data.BehavChannel;
    CaImgRawFN = data.CaImgRawFN;
    CaImgRawtime = data.CaImgRawtime;
    BehavRawFN = data.BehavRawFN;
    BehavRawtime = data.BehavRawtime;
    CaImgStartFN = data.CaImgStartFN;
    CaImgStartFtime = data.CaImgStartFtime;
    CaImgEndFN = data.CaImgEndFN;
    BehavStartFN = data.BehavStartFN;
    BehavStartFtime = data.BehavStartFtime;
    BehavEndFN = data.BehavEndFN;
    seg = data.seg;
    NeuStart = data.NeuStart;
    NeuEnd = data.NeuEnd;
    save(filename, "mouseID", "session", ...
        "CaImgChannel", "BehavChannel", ...
        "CaImgRawFN", "CaImgRawtime", ...
        "BehavRawFN", "BehavRawtime", ...
        "CaImgStartFN", "CaImgStartFtime", "CaImgEndFN", ...
        "BehavStartFN", "BehavStartFtime", "BehavEndFN", ...
        "seg", "NeuStart", "NeuEnd");
end

function defineInitData(data, filename)
    coordinates = data.coordinates;
    neuron = data.neuron;
    timestamp = data.timestamp;
    save(filename, "coordinates", "neuron", "timestamp");
end