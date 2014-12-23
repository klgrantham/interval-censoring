function [ymin, ymax] = extract_nth_ymin_ymax(n, data_dir)
% Extract ymin and ymax from an input data file

infile = fullfile(data_dir, sprintf('input_%06d.mat', n));

load(infile)

events = data.y;

if all(events(:, 2)==0 | events(:, 2)==1)  
    for i = 1:n
        tmp = events(i, 1);
        if events(i, 2)==1; events(i, 1)=tmp; events(i, 2)=tmp; end %fully observed
        if events(i, 2)==0; events(i, 1)=tmp; events(i, 2)=inf; end %right censored
    end
end

lo = events(:, 1);
ro = events(:, 2);
ally = unique([lo; ro(ro~=inf)]);
ymin = min(ally);
ymax = max(ally);

end