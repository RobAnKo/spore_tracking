function framovie = relabel_spores(framovie, Ms)
% This function needs: 
%   - framovie
%   - correspondence matrices
%
% and gives as output the same framovie, where framovie.cluster stores the relabeled spore masks

sporelabels = cellfun(@(x) label(x.spore, 1),framovie, 'UniformOutput', false);
origin = 1:max(sporelabels{1});
res_matrix = 1;

% produce correspondence matrix from frame 1 to last frame to determine complete tracks
end_mapping = 1;
for i = 1:length(Ms)
    end_mapping = end_mapping*Ms{i};
end

% complete tracks
remaining_tracks = unique(origin*end_mapping);

% discard all labels not corresponding to complete tracks
    % specifically for frame 1
    label_im = sporelabels{1};
    framovie{1}.cluster = label_im .* ismember(double(label_im), remaining_tracks);
for frame = 2:length(framovie)
    % second frame with unchanged labels
    label_im = double(sporelabels{frame});
    % matrix mapping IDs of first frame on the nth frame
    res_matrix = res_matrix * Ms{frame-1};
    % vectors of old labels where index is new label for value at index
    mapping = origin*res_matrix;
    out = zeros(size(label_im));
    for i = 1:size(label_im,1)
        for j = 1:size(label_im,2)
            before = label_im(i,j);
            if before ~= 0
                after = mapping(before);
                if ismember(after, remaining_tracks)
                    out(i,j) = after;
                end
            else
                continue
            end
        end
    end
    framovie{frame}.cluster = dip_image(out);
end
end    












