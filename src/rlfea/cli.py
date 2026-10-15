"""Command-line interface.

    rlfea audit GEOMETRY... [--unit mm]        geometry audit (no modification)
    rlfea prepare CONFIG [--runs DIR]          build mesh + solver input only (state: prepared)
    rlfea run CONFIG [--runs DIR]              prepare + solve + post-process + figures + report
    rlfea solve RUN_DIR [--restart DUMP --set etol=0.01]  (re)run or resume FEBio, then post-process
    rlfea post RUN_DIR                         derive fields, VTK export, checks, figures, report
    rlfea status RUN_DIR                       job state
    rlfea verify [--out DIR]                   analytical / benchmark verification suite
    rlfea capabilities                         FEBio availability and feature smoke tests
    rlfea gui [--runs DIR]                     interactive viewer (trame)
    rlfea export-web RUN_DIR --out f.html      standalone shareable HTML viewer of the solved fields
"""
from __future__ import annotations

import argparse
import json
import os
import sys


def _post(run_path):
    from . import runner
    from .report.latex import build_report
    runner.postprocess(run_path)
    runner.make_figures(run_path)
    out = build_report(run_path)
    print(json.dumps(out, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rlfea", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit"); a.add_argument("paths", nargs="+"); a.add_argument("--unit"); a.add_argument("--out")
    for name in ("prepare", "run"):
        p = sub.add_parser(name); p.add_argument("config"); p.add_argument("--runs", default="runs"); p.add_argument("--id")
        p.add_argument("--threads", type=int); p.add_argument("--timeout", type=float)
    s = sub.add_parser("solve"); s.add_argument("run_dir"); s.add_argument("--restart"); s.add_argument("--threads", type=int)
    s.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="restart override, e.g. etol=0.01 or max_retries=10 (recorded in the manifest)")
    sub.add_parser("post").add_argument("run_dir")
    sub.add_parser("status").add_argument("run_dir")
    v = sub.add_parser("verify"); v.add_argument("--out", default="verification_runs"); v.add_argument("--quick", action="store_true")
    sub.add_parser("capabilities")
    w = sub.add_parser("export-web"); w.add_argument("run_dir"); w.add_argument("--out", default="viewer.html")
    g = sub.add_parser("gui"); g.add_argument("--runs", default="runs"); g.add_argument("--port", type=int, default=8080)
    args = ap.parse_args(argv)

    if args.cmd == "audit":
        from .geometry.surface import audit_assembly, audit_body, load_body
        bodies = [b for p in args.paths for b in load_body(p, args.unit)]
        rep = {"bodies": {b.name: audit_body(b) for b in bodies}}
        if len(bodies) > 1:
            env = max((b for b in bodies if not b.empty), key=lambda b: abs(b.volume()))
            rep["assembly"] = audit_assembly(env, [b for b in bodies if b is not env])
            rep["assembly"]["envelope"] = env.name
        txt = json.dumps(rep, indent=1, default=str)
        if args.out:
            open(args.out, "w").write(txt)
        print(txt)
    elif args.cmd in ("prepare", "run"):
        from . import runner
        rd = runner.prepare(args.config, args.runs, args.id)
        print("prepared:", rd.path)
        if args.cmd == "run":
            res = runner.execute(rd, threads=args.threads, timeout=args.timeout)
            print("solver status:", res["status"], "last converged time:", res["convergence"].get("last_converged_time"))
            _post(rd.path)
            if res["status"] != "converged":
                sys.exit(2)
    elif args.cmd == "solve":
        from . import provenance as P, runner
        rd = P.RunDir(os.path.dirname(os.path.abspath(args.run_dir)), os.path.basename(os.path.abspath(args.run_dir)))
        ov = {}
        for kv in args.set:
            k, v = kv.split("=", 1)
            ov[k] = float(v) if any(c in v for c in ".e") else int(v)
        if ov and not args.restart:
            sys.exit("--set is only supported together with --restart")
        res = runner.execute(rd, threads=args.threads, restart=args.restart, overrides=ov)
        print(res["status"], "last converged time:", res["convergence"].get("last_converged_time"))
        _post(rd.path)
    elif args.cmd == "post":
        _post(args.run_dir)
    elif args.cmd == "status":
        from .provenance import read_status
        print(json.dumps(read_status(args.run_dir), indent=1))
    elif args.cmd == "verify":
        from .verification.suite import run_suite
        rep = run_suite(args.out, quick=args.quick)
        print(json.dumps({c["name"]: ("PASS" if c["pass"] else "FAIL") for c in rep["checks"]}, indent=1))
        sys.exit(0 if all(c["pass"] for c in rep["checks"]) else 1)
    elif args.cmd == "capabilities":
        from .solvers.capabilities import check
        print(json.dumps(check(), indent=1))
    elif args.cmd == "export-web":
        from .ui.web_export import export
        print(export(args.run_dir, args.out))
    elif args.cmd == "gui":
        from .ui.app import main as gui_main
        gui_main(args.runs, args.port)


if __name__ == "__main__":
    main()
