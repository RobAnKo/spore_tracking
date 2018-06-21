function fra = spore_segmentation(fra, cfg)

for i = 1:length(fra)
    imgRaw = fra{i}.bf/100;
    mask = imgRaw > cfg.INT_THRESH;
    mask = brmedgeobjs(mask);
    lab = label(mask, Inf, cfg.MIN_SIZE, cfg.MAX_SIZE);
    fra{i}.spore = lab > 0;
end