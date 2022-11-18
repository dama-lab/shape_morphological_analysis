function stat =  statsMetrics(array)
    stat.mean = mean(array);
    stat.std = std(array);
    stat.min = min(array);
    stat.max = max(array);
    stat.median = median(array);
end