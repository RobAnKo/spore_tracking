function dv_files = all_dv_to_position_lists(base_folder, positions_to_use, dvs_to_use)
%% read in all dv files in a folder and subdivide them by position

all_files = dir(fullfile(base_folder, '*lapse*.dv'));
split_names = arrayfun(@(x) split(x.name, {'_', 'lapse'}), all_files, 'UniformOutput', false);
for i = 1:size(split_names, 1)
    split_names{i}{2} = 'lapse';
end
all_dvs = arrayfun(@(x) str2double(x{1}{3}), split_names, 'UniformOutput', false);
all_positions = arrayfun(@(x) str2double(x{1}{4}), split_names, 'UniformOutput', false);
available_dvs = unique(cell2mat(all_dvs));
available_positions = unique(cell2mat(all_positions));

if string(positions_to_use) == "all"
    positions_to_use = available_positions;
else
    positions_to_use = str2num(positions_to_use)';
end% check whether dv-files for the chosen positions/frames exist in the folder
membership = ismember(positions_to_use, available_positions);
if any(~membership)
    missing = positions_to_use(~membership);
    error(['The positions ' num2str(missing) ' you choose are not given in the dv-files in the chosen folder'])
end

if string(dvs_to_use) == "all"
    dvs_to_use = available_dvs;
else
    dvs_to_use = str2num(dvs_to_use)';
end

membership = ismember(dvs_to_use, available_dvs);
if any(~membership)
    missing = dvs_to_use(~membership);
    error(['The dv-files ' num2str(missing) ' you choose are not given in the chosen folder'])
end


dv_files = cell(1,length(positions_to_use));
k = 1;
for position = positions_to_use'
    files = cell(1,length(dvs_to_use));
    j = 1;
    for file_index = 1:length(split_names)
        if str2double(split_names{file_index}{4}) == position && ismember(str2double(split_names{file_index}{3}), dvs_to_use)
            file = join(split_names{file_index}, '_');
            file = replace(file{1}, 'lapse_', 'lapse');
            files{j} = file;
            j = j+1;
        end
    end
    empty = cellfun(@(x) isempty(x), files);
    files = files(~empty);
    dv_files{k} = files;
    k = k + 1;
end
end
    


