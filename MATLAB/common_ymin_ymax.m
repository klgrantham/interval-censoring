function endpts = common_ymin_ymax(N, data_dir)

ymins = zeros(N,1);
ymaxs = zeros(N,1);
for i = 1:N
    [ymin, ymax] = extract_nth_ymin_ymax(i, data_dir);
    ymins(i) = ymin;
    ymaxs(i) = ymax;
end

max_ymins = max(ymins);
min_ymaxs = min(ymaxs);
endpts = [max_ymins; min_ymaxs];

end