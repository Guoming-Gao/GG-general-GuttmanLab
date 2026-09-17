"""Bounded, spawn-safe execution shared by notebook reconstruction pipelines."""
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from functools import wraps
from multiprocessing import get_context
import os

DEFAULT_MAX_WORKERS = 2


def worker_count(processing):
    value = processing.get("max_workers", DEFAULT_MAX_WORKERS)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("processing.max_workers must be a positive integer")
    return value


def initialize_worker():
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=1)


def pooled_batch(function):
    """Own one lazy process pool for an entire batch, including exceptional exits."""
    @wraps(function)
    def wrapped(config, *args, **kwargs):
        count = worker_count(config["processing"])
        if count == 1:
            return function(config, *args, executor=None, **kwargs)
        pool = ProcessPoolExecutor(
            max_workers=count, mp_context=get_context("spawn"), initializer=initialize_worker,
        )
        try:
            return function(config, *args, executor=pool, **kwargs)
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
    return wrapped


def execute_tasks(worker, tasks, *, executor, max_workers):
    """Yield completion order; cancel and settle a failed FOV before the next one."""
    if executor is None or max_workers == 1:
        for task in tasks:
            yield worker(task)
        return
    iterator = iter(tasks)
    pending = set()

    def fill():
        while len(pending) < 2 * max_workers:
            try:
                task = next(iterator)
            except StopIteration:
                break
            pending.add(executor.submit(worker, task))

    try:
        fill()
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                # Keep all other futures in pending until result() succeeds.
                result = future.result()
                pending.remove(future)
                yield result
            fill()
    finally:
        for future in pending:
            future.cancel()
        if pending:
            wait(pending)
