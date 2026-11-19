"""execute_with_resume: a run stopped by its time budget resumes from the last converged step."""
import json
import os

from rlfea import provenance as P
from rlfea import runner


def _rd(tmp_path, budget):
    rd = P.RunDir(str(tmp_path), "r")
    os.makedirs(rd.p("solver"), exist_ok=True)
    json.dump({"solver": {"timeout_s": budget}}, open(rd.p("config.resolved.json"), "w"))
    open(rd.p("solver", "model.dmp"), "w").write("dump")
    return rd


def test_resumes_after_interruption_until_converged(tmp_path, monkeypatch):
    rd = _rd(tmp_path, 18000)
    calls = []

    def fake_execute(rd_, threads=None, timeout=None, restart=None, overrides=None):
        calls.append((timeout, restart))
        return {"status": "interrupted" if len(calls) < 3 else "converged", "convergence": {}}

    monkeypatch.setattr(runner, "execute", fake_execute)
    res = runner.execute_with_resume(rd)
    assert res["status"] == "converged"
    assert calls[0][1] is None and calls[1][1].endswith("model.dmp") and calls[2][1].endswith("model.dmp")
    assert all(t <= 18000 for t, _ in calls)


def test_does_not_resume_a_failed_run(tmp_path, monkeypatch):
    rd = _rd(tmp_path, 18000)
    calls = []
    monkeypatch.setattr(runner, "execute", lambda *a, **k: calls.append(1) or {"status": "failed", "convergence": {}})
    assert runner.execute_with_resume(rd)["status"] == "failed"
    assert len(calls) == 1


def test_default_budget_is_five_hours():
    from rlfea.config import DEFAULTS
    assert DEFAULTS["solver"]["timeout_s"] == 18000
