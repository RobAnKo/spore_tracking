function framovie  = dvs2framovie(dvpath, dvfiles, UserInfo)
multiple_dvs = false;
% initialize frame and fquant cells
index = 0;
framovie{1} = frame();
fqmovie{1} = fquant();
% check: multiple files?
if iscell(dvfiles)
    multiple_dvs = true;
    len = length(dvfiles);
else
    len = 1;
end
% loop over all dv-files
for i = 1:len
    if multiple_dvs
        filename = fullfile(dvpath, dvfiles{i});
    else
        filename = fullfile(dvpath, dvfiles);
    end
    disp(['Currently processing ' filename]);
    %load the BF info from every dv-file
    dataDV = bfopen(filename);
    BFplanes = ReturnPlanes(dataDV, UserInfo.cfg.BF_CHANNEL, UserInfo);
    FLplanes = ReturnPlanes(dataDV, UserInfo.cfg.FL_CHANNELS, UserInfo);
    % check whether current dv-file contains multiple frames
    if length(BFplanes) > 1
        % loop over all frames in the dv-file
        for j = 1:length(BFplanes)
            index = index + 1;
            framovie{index} = frame();
            imgObject = dip_image(BFplanes{j}.ImageData);
            framovie{index}.bf = imgObject; %store BF signal in fra.bf
            
            framovie{index} = preprocess(framovie{index},UserInfo.cfg);
            framovie{index}.filename = sprintf('%s:%d', filename,j);
            
        end
    else
        index = index + 1;
        framovie{index} = frame();
        imgObject = dip_image(BFplanes{1}.ImageData);
        framovie{index}.bf = imgObject; %store BF signal in fra.bf
        framovie{index} = preprocess(framovie{index},UserInfo.cfg);
        framovie{index}.filename = sprintf('%s:%d', filename,1);
    end
end
end