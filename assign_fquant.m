function  fqarray = assign_fquant(dataDV,fra,BfPlanes,FlPlanes,UserInfo,quant_type,fqarray)

narginchk(6,7);

overwrite=UserInfo.Overwrite; 

assert(length(BfPlanes) == length(FlPlanes), 'QFTrack:internalError', 'BfPlanes length != FlPlanes length!');

if ~exist('fqarray', 'var')
    fqarray = cell(1,length(fra));
else
    assert(size(fqarray) == size(fra), 'QFTrack:internalError', 'fra length != fqarray length!');
end


for i=1:length(BfPlanes)
    FrameNr=BfPlanes{i}.FrameNr;
    
    if FrameNr <= length(fra)
    FrameNrsFl = structfun(@(x) x.FrameNr, FlPlanes{i});
    
    assert(all(FrameNr == FrameNrsFl));
    
    if ~overwrite
        if length(fqarray)>=FrameNr && ~isempty(fqarray{FrameNr})
            continue;
        end
    end
        
    ImgBfData=BfPlanes{i}.ImageData;
    
    if isa(ImgBfData,'cell')
        ImgBfData=ImgBfData{i};
    end
    ImgBfObject=dip_image(ImgBfData,class(ImgBfData));
    
    ImgFlObjects=struct();
    
    for channel_c = fieldnames(FlPlanes{i})'
        channel = char(channel_c);
        ImgFlData=FlPlanes{i}.(channel).ImageData;
%hack 6-30-14: put this part in, if you have fluorescent only pictures, take out again after the analysis is done
 %       if iscell(ImgFlData)
 %           ImgFlObjects.(channel) = dip_image(ImgFlData{1},class(ImgFlData{1}));
 %       else
 %           ImgFlObjects.(channel) = dip_image(ImgFlData,class(ImgFlData));
 %       end
        %this is the standard code for processing images with more than
        %one channel
           ImgFlObjects.(channel) = dip_image(ImgFlData,class(ImgFlData));
    end
    

        if ~isempty(fra{FrameNr})
          %  fq = measurecells(fquant(),fra{FrameNr},ImgBfObject, ImgFlObjects, dataDV{2});
          %fq = measure_all(fquant(),fra{FrameNr},ImgBfObject, ImgFlObjects, dataDV{2});
          fq = generate_fquant(fquant(),fra{FrameNr}, ImgFlObjects, dataDV{2},quant_type);
    
        else
            fq=[];
        end
    %Ilka 21.5.2012
    %maybe this is causing trouble, but otherwise we run at risk to have an
    %incomplete time vector, leave in for the time being, otherwise NaN
    %problem?
        if isempty(fq)
            continue
        end
    
        bfplanenr=BfPlanes{i}.PlaneNr;
        fq.general=SetLogData(dataDV,bfplanenr,fq.general,UserInfo);
    
        fqarray{FrameNr} = fq;
    end
    
    
end

disp('Assign_fquant: fq done');
