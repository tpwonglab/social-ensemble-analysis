%% Main File to simulate entire Social Valence data pipeline
main()

function main
    setupSource()
    createTestData()
end

function setupSource
    disp("Run CNMF script to view Sources2D data");
    run CNMF_E/cnmfe_setup.m;
end

function createTestData
    disp("Create simulation configuration file");
    data = load("CSDS/1035_Def8.mat");
    defineInitExp(data.segment, "init_two_seg.mat");
    defineInitData(data, "init_two_data.mat");
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
    binwidth = data.binwidth;
    save("data/" + filename, "mouseID", "session", ...
        "CaImgChannel", "BehavChannel", ...
        "CaImgRawFN", "CaImgRawtime", ...
        "BehavRawFN", "BehavRawtime", ...
        "CaImgStartFN", "CaImgStartFtime", "CaImgEndFN", ...
        "BehavStartFN", "BehavStartFtime", "BehavEndFN", ...
        "seg", "NeuStart", "NeuEnd", "binwidth");
end

function defineInitData(data, filename)
    coordinates = data.coordinates;
    neuron = data.neuron;
    timestamp = data.timestamp;
    save("data/" + filename, "coordinates", "neuron", "timestamp");
end