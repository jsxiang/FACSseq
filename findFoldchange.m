function [foldchange,switchhits]=findFoldchange(mulig,munolig,switchthreshold)
    foldchange=10.^(mulig-munolig);
    switchhits=foldchange>switchthreshold;
end