"""
[COPY with two additions for the Figure 6 v3 add-on, 2026-10-03: retry on transient server errors and explicit
memory release per variant. Prediction settings are unchanged from the original alphagenome_perbase_MRAS.py.]

AlphaGenome PER-BASE (base-resolution) re-run for the MRAS locus.

Unlike score_variant_api.py (which returns *aggregated* CenterMask/GeneMask scores that
cannot separate the four MRAS promoters — their TSS are within ~1.5 kb), this script keeps
the base-resolution REF and ALT track profiles so the effect can be read at EACH promoter TSS.

For every variant it calls predict_variant() over a 1 Mb interval, then slices every requested
modality's REF/ALT TrackData down to a small window covering all four MRAS promoters and saves:
  per_snp/<snp_id>.npz  ->  arrays  ref_<MOD>, alt_<MOD>  (shape: n_positions x n_tracks),
                            plus  positions (genomic, 0-based) and window start/end.
  meta/<MOD>_tracks.parquet  ->  track metadata for that modality (written once).

Usage:
  python alphagenome_perbase_MRAS.py <snplist.txt> <out_dir> [MOD1,MOD2,...]
    snplist : one snp_id per line (chr3_138348784_C_T_b38) — or TSV with a 'snp' column
    out_dir : output directory (per_snp/ and meta/ created inside)
    MODs    : comma list of modalities (default: CAGE,PROCAP,DNASE,ATAC)
              valid: CAGE PROCAP DNASE ATAC CHIP_HISTONE CHIP_TF RNA_SEQ SPLICE_SITES

MRAS promoters (GRCh38, +):  P45111 TSS 138347648 | P45112 138348579 | P45113 138348850 | P45114 138349146
Window default 138345000-138352000 contains all four plus flanks.
"""
import os, sys, time, gc
import numpy as np
import pandas as pd

from alphagenome.models import dna_client
from alphagenome.data import genome

API_KEY   = os.environ['ALPHAGENOME_API_KEY']   # set in the environment; the key is not stored in the repository
SUPPORTED = sorted(dna_client.SUPPORTED_SEQUENCE_LENGTHS.values())  # 16KB,100KB,500KB,1MB
WIN_CHR   = 'chr3'
WIN_START = 138_345_000      # 0-based inclusive start of the saved per-base window
WIN_END   = 138_352_000      # exclusive end  (7000 bp -> covers all 4 MRAS promoters + flanks)

def parse_snps(fp):
    if fp.endswith(('.tsv', '.csv')):
        sep = '\t' if fp.endswith('.tsv') else ','
        df = pd.read_csv(fp, sep=sep)
        return df['snp'].astype(str).unique().tolist()
    with open(fp) as fh:
        return [l.strip() for l in fh if l.strip()]

def main():
    snplist = sys.argv[1]
    out_dir = sys.argv[2]
    mods = (sys.argv[3].split(',') if len(sys.argv) > 3
            else ['CAGE', 'PROCAP', 'DNASE', 'ATAC'])
    out_types = [getattr(dna_client.OutputType, m) for m in mods]

    per_snp = os.path.join(out_dir, 'per_snp'); os.makedirs(per_snp, exist_ok=True)
    meta_d  = os.path.join(out_dir, 'meta');    os.makedirs(meta_d,  exist_ok=True)

    snps = parse_snps(snplist)
    done = {f[:-4] for f in os.listdir(per_snp) if f.endswith('.npz')}
    todo = [s for s in snps if s not in done]
    print(f'modalities={mods}  window={WIN_CHR}:{WIN_START}-{WIN_END} ({WIN_END-WIN_START} bp)')
    print(f'total={len(snps)}  done={len(done)}  todo={len(todo)}')

    print('connecting to AlphaGenome API...')
    model = dna_client.create(API_KEY)
    print('connected.')

    def predict_with_backoff(interval, variant, tries=8):
        """Retry on API quota (RESOURCE_EXHAUSTED): the AlphaGenome MB/minute limit resets
        each minute, so back off ~that long and retry instead of dropping the variant."""
        delay = 20
        for attempt in range(1, tries + 1):
            try:
                return model.predict_variant(interval=interval, variant=variant,
                                             requested_outputs=out_types, ontology_terms=None)
            except Exception as e:
                msg = str(e)
                if 'RESOURCE_EXHAUSTED' in msg and attempt < tries:
                    print(f'      quota hit (attempt {attempt}), sleeping {delay}s...', flush=True)
                    time.sleep(delay); delay = min(delay + 15, 75)
                elif any(k in msg for k in ('StatusCode.INTERNAL', 'StatusCode.UNAVAILABLE', 'DEADLINE_EXCEEDED')) and attempt < tries:
                    # ADDED (Figure 6, 2026-10-03): transient server-side errors on 1 Mb requests -> back off and retry
                    print(f'      server error (attempt {attempt}): {msg.splitlines()[0][:80]} ... sleeping {delay * 3}s', flush=True)
                    gc.collect(); time.sleep(delay * 3); delay = min(delay + 15, 75)
                else:
                    raise

    fixed_seqlen = int(os.environ.get('AG_SEQLEN', '0'))  # e.g. 524288 = uniform 500KB
    def choose_interval(variant):
        """Interval that contains BOTH the variant and the promoter window, centred on their union.
        If AG_SEQLEN is set (e.g. 524288 = 500KB), use that fixed length for ALL variants
        (uniform => cross-variant comparable, matches the puQTL promoter+/-500kb cis-window);
        it is only bumped up if a variant+window span exceeds it. Otherwise adaptive minimal."""
        vi = variant.reference_interval
        lo = min(vi.start, WIN_START); hi = max(vi.end, WIN_END)
        need = hi - lo
        if fixed_seqlen and fixed_seqlen >= need:
            W = fixed_seqlen
        else:
            W = next((w for w in SUPPORTED if w >= need), SUPPORTED[-1])
        return genome.Interval(vi.chromosome, lo, hi).resize(W), W

    meta_written = set()
    for i, snp in enumerate(todo):
        t0 = time.time()
        try:
            variant  = genome.Variant.from_str(snp, variant_format=genome.VariantFormat.GTEX)
            interval, W = choose_interval(variant)
            vo = predict_with_backoff(interval, variant)
            save = {}
            for m in mods:
                key = m.lower()
                ref_td = getattr(vo.reference, key)
                alt_td = getattr(vo.alternate, key)
                iv = ref_td.interval
                # base index of window inside the returned interval
                lo = WIN_START - iv.start
                hi = WIN_END   - iv.start
                lo = max(lo, 0); hi = min(hi, ref_td.values.shape[0])
                save[f'ref_{m}'] = ref_td.values[lo:hi].astype(np.float32)
                save[f'alt_{m}'] = alt_td.values[lo:hi].astype(np.float32)
                if f'pos_{m}' not in save:
                    save[f'pos_{m}'] = np.arange(iv.start + lo, iv.start + hi, dtype=np.int64)
                # track metadata (written once per modality; skip if present -> safe for parallel jobs)
                mp = os.path.join(meta_d, f'{m}_tracks.parquet')
                if m not in meta_written and not os.path.exists(mp):
                    ref_td.metadata.copy().to_parquet(mp)
                    meta_written.add(m)
            save['win_start'] = np.int64(WIN_START)
            save['win_end']   = np.int64(WIN_END)
            save['seqlen']    = np.int64(W)
            np.savez_compressed(os.path.join(per_snp, f'{snp}.npz'), **save)
            print(f'  [{i+1}/{len(todo)}] {snp}  seqlen={W}  ({time.time()-t0:.1f}s)', flush=True)
        except Exception as e:
            print(f'  [{i+1}/{len(todo)}] {snp}  ERROR ({time.time()-t0:.1f}s): {str(e).splitlines()[0][:200]}', flush=True)
        finally:
            # ADDED (Figure 6, 2026-10-03): free the multi-GB prediction before the next variant
            vo = None; save = None; gc.collect()

    print('done. per-base npz in', per_snp)

if __name__ == '__main__':
    main()
