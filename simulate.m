clear main;
clc;
main()

function main
    addpath("functions");    

    % setupSource()
    % 
    % segmentFilename = "data/init_two_seg.mat";
    % dataFilename = "data/init_two_data.mat";
    % isSingleMouse = 0;
    % createTestData(segmentFilename, dataFilename);
    % 
    % outputFilename = binning(segmentFilename, dataFilename, isSingleMouse);
    % 
    % csi(outputFilename);

    ensemble(outputFilename, 1);
end

function setupSource
    run CNMF_E/cnmfe_setup.m;
end

function createTestData(segmentFilename, dataFilename)
    if exist(segmentFilename, "file") && exist(dataFilename, "file")
        return
    end
    data = load("CSDS/1035_Def8.mat");
    defineInitExp(data.segment, segmentFilename);
    defineInitData(data, dataFilename);
end

function defineInitExp(data, filename)
    mouseID = "1";
    session = "test";
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
    save("data/" + filename, "mouseID", "session", ...
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