# Fetch current Starlink public-ephemeris vintage, compute cross-vintage covariance
# realism statistics, persist compact artifacts. Runs on GitHub Actions (ubuntu) or locally.
# Source: https://api.starlink.com/public-files/ephemerides/ (public, no auth).
# Raw files (~2 MB each, ~7000 sats ~= 14 GB per vintage) are processed in flight and
# discarded; only compact derived products are kept (see README).
import gzip, json, os, re, time, urllib.request, datetime
from concurrent.futures import ThreadPoolExecutor
import numpy as np

BASE = "https://api.starlink.com/public-files/ephemerides/"
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Chrome/131.0 satellite-ephemeris-archive/1.0"}
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUBSET_N = int(os.environ.get("SUBSET_N", "200"))   # sats with detailed sampled states
SUBSET_STEP = 30                                    # index step (=minutes; file cadence 60 s)
REF_OFFS_H = np.array([2, 6, 12, 24, 36, 48])       # sparse ref epochs, all sats
CONC = int(os.environ.get("CONC", "8"))
MANEUVER_M = 1000.0
BUCKET_EDGES_H = np.array([0, 6, 12, 24, 36, 48])

def now_utc(): return datetime.datetime.now(datetime.timezone.utc)

def get(url, timeout=180, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))

EPOCH_RE = re.compile(r"^\d{13}\.\d")
def parse_meme(text):
    """Modified ITC text -> (t_unix [n], rv [n,6] km|km/s, cov [n,6,6] symmetric)."""
    lines = text.splitlines()
    t = np.empty(len(lines) // 4 + 2); rv = np.empty((len(t), 6)); cov = np.empty((len(t), 6, 6))
    k = 0
    for i in range(len(lines)):
        toks = lines[i].split()
        if toks and EPOCH_RE.match(toks[0]) and len(toks) == 7:
            ep = toks[0]
            doy = datetime.date(int(ep[:4]), 1, 1) + datetime.timedelta(days=int(ep[4:7]) - 1)
            base = datetime.datetime(doy.year, doy.month, doy.day, int(ep[7:9]),
                                     int(ep[9:11]), int(ep[11:13]),
                                     tzinfo=datetime.timezone.utc)
            t[k] = base.timestamp() + float(ep[13:])
            rv[k] = [float(x) for x in toks[1:7]]
            c = [float(x) for l2 in lines[i + 1:i + 4] for x in l2.split()]
            ii, jj = np.triu_indices(6)
            cov[k][ii, jj] = cov[k][jj, ii] = c
            k += 1
    return t[:k], rv[:k], cov[:k]

def hermite_pos(t, rv, tq):
    j = np.clip(np.searchsorted(t, tq) - 1, 0, len(t) - 2)
    h = t[j + 1] - t[j]; s = ((tq - t[j]) / h)[:, None]      # (n,1) broadcasts with (n,3)
    p0, v0, p1, v1 = rv[j, :3], rv[j, 3:], rv[j + 1, :3], rv[j + 1, 3:]
    return ((2*s**3 - 3*s**2 + 1) * p0 + (s**3 - 2*s**2 + s) * (h[:, None] * v0)
            + (-2*s**3 + 3*s**2) * p1 + (s**3 - s**2) * (h[:, None] * v1))

def main():
    ts = now_utc().strftime("%Y%m%dT%H%M%SZ")
    rep = {"run": ts}
    for d in ("manifests", "pairs", "states_subset/" + ts[:6], "state"):
        os.makedirs(f"{ROOT}/{d}", exist_ok=True)

    manifest = get(BASE + "MANIFEST.txt").decode()
    with gzip.open(f"{ROOT}/manifests/{ts}_MANIFEST.txt.gz", "wt") as f:
        f.write(manifest)
    files = [f for f in manifest.split() if f.strip()]
    rep["manifest_files"] = len(files)
    sat_of = lambda fn: re.match(r"MEME_\d+_(STARLINK-[\dA-Z]+)_", fn).group(1)
    sat2file = {sat_of(fn): fn for fn in files}
    sats_all = sorted(sat2file)
    if os.environ.get("LIMIT_SATS"):          # smoke-test mode: first N sats only
        sats_all = sats_all[:int(os.environ["LIMIT_SATS"])]
    subset = set(sats_all[::max(1, len(sats_all) // SUBSET_N)])
    rep["sats_total"], rep["sats_subset"] = len(sats_all), len(subset)

    prev = None
    if os.path.exists(f"{ROOT}/state/refs.npz"):
        z = np.load(f"{ROOT}/state/refs.npz")
        prev = {k: z[k] for k in z.files}
        prev_index = {s: i for i, s in enumerate(prev["sat_ids"].tolist())}
        do_pairs = prev["vintage"][0] != files[0]
    else:
        prev, prev_index, do_pairs = None, {}, False
    rep["do_pairs"] = do_pairs

    def handle(sat):
        txt = get(BASE + sat2file[sat]).decode()
        t, rv, cov = parse_meme(txt)
        if len(t) < 10:
            return None
        idx = np.clip(np.searchsorted(t, t[0] + REF_OFFS_H * 3600.0), 0, len(t) - 1)
        out = {"sat": sat, "ref_t": t[idx], "ref_rv": rv[idx].astype(np.float32),
               "ref_cov": cov[idx].astype(np.float32), "t0": t[0]}
        if sat in subset:
            out["sub"] = (t[::SUBSET_STEP], rv[::SUBSET_STEP].astype(np.float32),
                          cov[::SUBSET_STEP].astype(np.float32))
        if do_pairs and sat in prev_index:
            si = prev_index[sat]
            pt, prv, pcov = prev["ref_t"][si], prev["ref_rv"][si], prev["ref_cov"][si]
            pt0 = prev["t0"][si]
            ok = np.isfinite(pt)
            tq = pt[ok]
            m = (tq > t[0] + 2700) & (tq < t[-1] - 2700)
            tq = tq[m]
            if len(tq):
                p2 = hermite_pos(t, rv, tq)
                dr_km = p2 - prv[ok][m][:, :3]
                dr_m = np.linalg.norm(dr_km, axis=1) * 1e3
                P_m = pcov[ok][m][:, :3, :3] * 1e6
                chi2 = np.einsum("ni,nij,nj->n", dr_km, np.linalg.inv(P_m), dr_km) * 1e6
                age = tq - pt0
                b = []
                for lo, hi in zip(BUCKET_EDGES_H[:-1], BUCKET_EDGES_H[1:]):
                    sel = (age > lo * 3600) & (age <= hi * 3600)
                    b.append(float(np.median(dr_m[sel])) if sel.any() else np.nan)
                out["pair"] = (sat, len(tq), float(np.median(dr_m)), float(np.median(chi2)),
                               *b, 1 if np.median(dr_m) > MANEUVER_M else 0)
        return out

    t00 = time.time()
    ref_rows, sub_rows, pair_rows, n_fail = [], [], [], 0
    with ThreadPoolExecutor(max_workers=CONC) as ex:
        for r in ex.map(handle, sats_all):
            if r is None:
                n_fail += 1
                continue
            ref_rows.append(r)
            if "sub" in r:
                sub_rows.append((r["sat"],) + r["sub"])
            if "pair" in r:
                pair_rows.append(r["pair"])
    rep.update(parse_seconds=round(time.time() - t00, 1), parsed=len(ref_rows), failed=n_fail)
    if n_fail > len(sats_all) * 0.1:
        rep["WARN"] = "more than 10% of files failed"

    # persist refs for next run
    S, R = len(ref_rows), len(REF_OFFS_H)
    sat_arr = np.array([r["sat"] for r in ref_rows])
    ref_t = np.full((S, R), np.nan); ref_rv = np.zeros((S, R, 6), np.float32)
    ref_cov = np.zeros((S, R, 6, 6), np.float32); t0 = np.zeros(S)
    for i, r in enumerate(ref_rows):
        ref_t[i], ref_rv[i], ref_cov[i], t0[i] = r["ref_t"], r["ref_rv"], r["ref_cov"], r["t0"]
    np.savez(f"{ROOT}/state/refs.npz", sat_ids=sat_arr, ref_t=ref_t, ref_rv=ref_rv,
             ref_cov=ref_cov, t0=t0, vintage=np.array([files[0]]))
    open(f"{ROOT}/state/last_manifest.txt", "w").write(manifest)

    if sub_rows:
        n_max = max(len(x[1]) for x in sub_rows)
        def pad(arrs):
            out = np.full((len(arrs), n_max) + np.shape(arrs[0])[1:], np.nan, arrs[0].dtype)
            for i, a in enumerate(arrs):
                out[i, :len(a)] = a
            return out
        np.savez(f"{ROOT}/states_subset/{ts[:6]}/{ts}.npz",
                 sat_ids=np.array([s for s, *_ in sub_rows]),
                 t=pad([x[1] for x in sub_rows]),
                 rv=pad([x[2] for x in sub_rows]),
                 cov=pad([x[3] for x in sub_rows]))
    rep["subset_saved"] = len(sub_rows)

    if pair_rows:
        with gzip.open(f"{ROOT}/pairs/{ts}.csv.gz", "wt") as f:
            f.write("sat,n,med_dr_m,med_chi2,b_0_6h,b_6_12h,b_12_24h,b_24_36h,b_36_48h,flag\n")
            for row in pair_rows:
                f.write(",".join(repr(x) for x in row) + "\n")
        drs = sorted(r[2] for r in pair_rows)
        nm = [r[2] for r in pair_rows if r[2] <= MANEUVER_M]
        rep.update(pair_sats=len(pair_rows),
                   med_dr_nonmaneuver_m=round(float(np.median(nm)), 1) if nm else None,
                   maneuver_like=sum(1 for r in pair_rows if r[9] == 1))
        rep["dr_p50_p90_p99_m"] = [round(drs[len(drs) // 2], 1),
                                   round(drs[int(len(drs) * .9)], 1),
                                   round(drs[int(len(drs) * .99)], 1)]

    sizes = {d: round(sum(os.path.getsize(os.path.join(dp, fn)) for dp, _, fns in
                          os.walk(f"{ROOT}/{d}") for fn in fns) / 1e6, 1)
             for d in ("manifests", "pairs", "states_subset", "state")}
    rep["sizes_mb"] = sizes
    with open(f"{ROOT}/runs.log", "a") as f:
        f.write(json.dumps(rep) + "\n")
    print(json.dumps(rep, indent=1))

if __name__ == "__main__":
    main()
