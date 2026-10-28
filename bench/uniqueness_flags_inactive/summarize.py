import csv, statistics as st, sys
for path in sys.argv[1:]:
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    print(path.split('/')[-1])
    med = {}
    for v in ('base', 'flags'):
        r = [x for x in rows if x['binary'] == v]
        w = [float(x['wall_s']) for x in r]; u = [float(x['user_s']) for x in r]
        med[v] = (st.median(w), st.median(u))
        print(f"  {v:5s} n={len(r)} wall median {st.median(w):7.2f}s (min {min(w):.2f}, max {max(w):.2f})  "
              f"user median {st.median(u):7.2f}s (min {min(u):.2f}, max {max(u):.2f})  rss {max(int(x['maxrss_kb']) for x in r)/1e6:.2f} GB")
    print(f"  flags vs base: wall {100*(med['flags'][0]/med['base'][0]-1):+.2f}%  user {100*(med['flags'][1]/med['base'][1]-1):+.2f}%")
