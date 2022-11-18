function [surface_area,csv_cell] = calculate_surface_area(full_volume_path,thickness_path,surface_area_path)
        % full volume
        [volume, csv_cell] = extract_feature_from_csv(full_volume_path);
        % full thickness
        thickness = extract_feature_from_csv(thickness_path);
        % calculate the surface area
        surface_area = volume ./ thickness;
        % convert matrix to cell array
        csv_cell(3:end,3:end) = num2cell(surface_area);
        % save_surface_area_file
        if exist('surface_area_path','var')
            writetable(cell2table(csv_cell),surface_area_path,'WriteVariableNames',false);
        end
        
end

function [feature,csv_cell] = extract_feature_from_csv(csv_path)
    csv_cell = readCsvToCell(csv_path);
    feature = cellfun(@str2double,csv_cell(3:end,3:end));
end
