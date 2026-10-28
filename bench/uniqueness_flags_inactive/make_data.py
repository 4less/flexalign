"""Reference (first genome of the 60 largest Bacteroides species) and simulated 150 bp pairs."""
import gzip, random, sys
from collections import OrderedDict
mapfile, out, npairs = sys.argv[1], sys.argv[2], int(sys.argv[3])
random.seed(7)
species = OrderedDict()
for line in open(mapfile):
    acc, path, tax = line.rstrip('\n').split('\t')
    species.setdefault(tax.split(';')[-1], []).append(path)
picked = sorted(species.items(), key=lambda kv: -len(kv[1]))[:60]
seqs = []
with open(f'{out}/reference.fna', 'w') as ref, open(f'{out}/sequence_species.tsv', 'w') as tsv:
    for sp, paths in picked:
        name, buf = None, []
        for line in gzip.open(paths[0], 'rt'):
            if line.startswith('>'):
                if name: seqs.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0], []
                ref.write(line.split()[0] + '\n'); tsv.write(f'{name}\t{sp}\n')
            else:
                buf.append(line.strip().upper()); ref.write(line)
        seqs.append((name, ''.join(buf)))
seqs = [s for s in seqs if len(s[1]) > 2000]
weights = [len(s) for _, s in seqs]
comp = str.maketrans('ACGT', 'TGCA')
def mut(r):
    r = list(r)
    for i in range(len(r)):
        if random.random() < 0.01: r[i] = random.choice('ACGT')
    return ''.join(r)
q = 'I' * 150
with open(f'{out}/reads_R1.fq', 'w') as f1, open(f'{out}/reads_R2.fq', 'w') as f2:
    for i, (name, s) in enumerate(random.choices(seqs, weights, k=npairs)):
        frag = int(random.gauss(400, 50)); frag = max(200, min(frag, len(s) - 1))
        p = random.randrange(0, len(s) - frag)
        a, b = s[p:p + 150], s[p + frag - 150:p + frag].translate(comp)[::-1]
        if random.random() < 0.5: a, b = b, a
        f1.write(f'@r{i}/1\n{mut(a)}\n+\n{q}\n'); f2.write(f'@r{i}/2\n{mut(b)}\n+\n{q}\n')
