def TCR_EM_counts(unique_counts, counts_old, t_indices, thresh, max_iters):
    working = True
    iters = 0
    while working:
        iters += 1
        counts = list(unique_counts)
        print('update ' + str(iters))
        #print('first element of counts: ' + str(counts[0]))
        #print('first element of counts_old: ' + str(counts_old[0]))
        for idx in t_indices:
            vals = [counts_old[x-1] for x in idx]
            s = sum(vals)
            for ii in range(len(idx)):
                counts[idx[ii]-1] += vals[ii]/s
        diff = max([abs(counts[ii]-counts_old[ii]) for ii in range(len(counts))])
        #print(diff)
        if (diff < thresh) or (iters >= max_iters):
            working = False
        else:
            counts_old = list(counts)
    return counts

