clear main;
clc;
main()

function main
    addpath("functions");
    setupSource()
    segmentFilename = "data/init_one_seg.mat";
    dataFilename = "data/init_one_data.mat";
    isSingleMouse = 1;
    createTestData(segmentFilename, dataFilename);
    binning(segmentFilename, dataFilename, isSingleMouse);
end

function setupSource
    disp("Run CNMF script to view Sources2D data");
    run CNMF_E/cnmfe_setup.m;
end

function createTestData(segmentFilename, dataFilename)
    if exist(segmentFilename, "file") && exist(dataFilename, "file")
        return
    end
    disp("Create simulation configuration file");
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