from concurrent.futures import Future, ThreadPoolExecutor
from io import StringIO
from threading import Event
from unittest.mock import patch

import pytest
from rich.console import Console

from sacdpy.execution import execute_tasks, pooled_batch, worker_count
from sacdpy.progress import PipelineProgress, progress_items, report_progress, stage


def test_worker_count_validation_and_default():
    assert worker_count({}) == 2
    assert worker_count({"max_workers": 1}) == 1
    for value in (0, -1, True, 1.5, "3"):
        with pytest.raises(ValueError):
            worker_count({"max_workers": value})


def test_out_of_order_completion_and_fewer_tasks_than_workers():
    released = Event()
    def worker(i):
        if i == 0:
            assert released.wait(5)
        return i
    with ThreadPoolExecutor(3) as pool:
        results = execute_tasks(worker, [0, 1], executor=pool, max_workers=3)
        assert next(results) == 1
        released.set()
        assert list(results) == [0]


@pytest.mark.parametrize("error", [ValueError("worker failed"), KeyboardInterrupt()])
def test_failure_cancels_pending_and_settles_running_jobs(error):
    class Pool:
        def __init__(self):
            self.futures = []
        def submit(self, worker, task):
            future = Future()
            if not self.futures:
                future.set_exception(error)
            self.futures.append(future)
            return future
    pool = Pool()
    # Pending futures need the executor's cancelled notification for wait().
    with patch("sacdpy.execution.wait", side_effect=lambda fs, **kw: ({f for f in fs if f.done()}, {f for f in fs if not f.done()})) as settled:
        with pytest.raises(type(error)):
            list(execute_tasks(lambda x: x, range(20), executor=pool, max_workers=3))
    assert len(pool.futures) == 6
    assert all(f.done() for f in pool.futures)
    assert settled.call_count == 2


def test_pool_shutdown_even_when_batch_raises():
    @pooled_batch
    def failing(config, *, executor):
        raise ValueError("preflight failure")
    with patch("sacdpy.execution.ProcessPoolExecutor") as factory:
        with pytest.raises(ValueError):
            failing({"processing": {}})
        factory.return_value.shutdown.assert_called_once_with(wait=True, cancel_futures=True)


def test_progress_resumed_failed_and_dataset_identity():
    console = Console(file=StringIO(), force_terminal=False)
    with PipelineProgress(9, console=console, auto_refresh=False) as progress:
        progress({"event": "reconstruction_done", "relative_fov": "FOV", "dataset_name": "A", "completed_in_fov": 1, "total_in_fov": 3})
        progress({"status": "failed", "relative_fov": "FOV", "dataset_name": "A", "movie_count": 3, "error": "worker"})
        progress({"status": "resumed_manifest", "relative_fov": "FOV", "dataset_name": "B", "movie_count": 3})
        progress({"status": "resumed_manifest", "relative_fov": "other", "movie_count": 3})
        assert progress.progress.tasks[progress.overall].completed == 7
        assert progress.resumed == 2
        assert progress.failed == 1
        assert progress.failed_work == 2


def test_all_resumed_and_interrupted_progress():
    progress = PipelineProgress(6, console=Console(file=StringIO()), auto_refresh=False)
    with pytest.raises(KeyboardInterrupt):
        with progress:
            for fov in ("a", "b"):
                progress({"status": "resumed_manifest", "relative_fov": fov, "reconstruction_count": 3})
            assert progress.progress.tasks[progress.overall].completed == 6
            raise KeyboardInterrupt()
    assert progress.cancelled


def test_serial_stage_callbacks_inherit_and_reset():
    events = []
    @report_progress
    def inner():
        for _ in progress_items([1, 2], "nuclei"):
            pass
    @report_progress
    def outer():
        stage("load model")
        inner()
    outer(progress_callback=events.append)
    assert events[0]["total"] is None
    assert events[-1]["completed"] == events[-1]["total"] == 2
    stage("outside")
    assert len(events) == 4
