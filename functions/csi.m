function csi(filename)
    if isempty(filename)
        quit(2);
    end
    data = load(filename);
    disp(data);
end