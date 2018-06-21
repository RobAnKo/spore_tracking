function [fq, framovie] = dv2fq(framovie, dvpath, UserInfo)
%% This function imports a dv-file and quantifies its content based on the arguments given in UserInfo

quant_type = UserInfo.cfg.quant_type;




disp(['Currently processing ' dvpath]);

%% import data from every dv-file
dataDV = bfopen(dvpath);





% extract BF info 
BFplanes = ReturnPlanes(dataDV, UserInfo.cfg.BF_CHANNEL, UserInfo);
% extract FL info
for n = 1:length(UserInfo.cfg.FL_CHANNELS) %for all channels asked
    % get channel number
    Fl_channel = UserInfo.cfg.FL_CHANNELS(n);
    % get the image information of the Flchannel
    cur_flplanes = ReturnPlanes(dataDV,Fl_channel,UserInfo);
    
    if isempty(cur_flplanes)
        error('QFTrack:internalError', 'no fl planes returned');
    end
    
    % extract the name of fluorescent agent
    PlaneFl = num2str(cur_flplanes{1}.PlaneNr);
    stringImage = ['Image ' PlaneFl '. EX filter'];
    flchannel = dataDV{2}.get(stringImage);
    if isempty(flchannel)
        flchannel = 'unknown';
    else
        flchannel = strrep(flchannel, '-', '_'); %varname does not accept '-'s
    end
            
    for i = 1:length(cur_flplanes)
        FlPlanes{i}.(flchannel) = cur_flplanes{i};
    end
    
end
%% transfer data into fra object
fra = cell(1,length(BFplanes));
for i=1:length(BFplanes)
    FrameNr = BFplanes{i}.FrameNr;
    ImgData = BFplanes{i}.ImageData;
    ImgFlData = structfun(@(x) x.ImageData, FlPlanes{i}, 'UniformOutput', false);
        
    if isa(ImgData,'cell')
        ImgData = ImgData{1};
    end
    
    ImgObject=dip_image(ImgData,class(ImgData));
    ImgFlObject = structfun(@dip_image, ImgFlData, 'UniformOutput', false);
    
    disp(['Frame ' num2str(FrameNr)]);
    
    % frame initialization
    fra{FrameNr} = frame();
    fra{FrameNr}.bf = ImgObject;
    % preprocessing
    fra{FrameNr} = preprocess(fra{FrameNr},UserInfo.cfg);
    fra{FrameNr}.fl = ImgFlObject;
    fra{FrameNr}.filename = sprintf('%s:%d', dvpath,FrameNr);
    fra{FrameNr}.cfg = UserInfo.cfg;
    currName = [dvpath '_' num2str(FrameNr, '%02d\n')];
    
end
% segmentation
fra = spore_segmentation(fra,UserInfo.cfg);
% quantification
framovie = [framovie fra];
fq = assign_fquant(dataDV, fra, BFplanes, FlPlanes, UserInfo, quant_type); 
end