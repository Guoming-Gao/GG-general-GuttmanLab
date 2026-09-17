"""Parent-process Rich rendering and optional callbacks for serial analysis stages."""
from contextvars import ContextVar
from functools import wraps

from rich.progress import (
    BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn,
    TimeElapsedColumn, TimeRemainingColumn,
)


def make_progress(**kwargs):
    return Progress(
        SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
        BarColumn(), MofNCompleteColumn(), TimeElapsedColumn(), TimeRemainingColumn(),
        **kwargs,
    )


class FOVEvents:
    """Track failure stage independently of whether a UI callback is supplied."""
    def __init__(self, callback):
        self.callback = callback
        self.failure_stage = "reconstruction"

    def __call__(self, event):
        if event.get("event") == "validation_started":
            self.failure_stage = "output_validation_or_publication"
        if self.callback is not None:
            self.callback(event)


class PipelineProgress:
    """Counts computed/resumed work without marking failed work as successful."""
    def __init__(self, total=None, *, description="All reconstructions", **kwargs):
        self.progress = make_progress(**kwargs)
        self.overall = self.progress.add_task(description, total=total)
        self.current = self.progress.add_task("Waiting", total=None, visible=False)
        self.counts = {}
        self.resumed = 0
        self.failed = 0
        self.failed_work = 0
        self.cancelled = False
        self.description = description
        self.current_key = None
        self.stage_names = {}

    def __enter__(self):
        self.progress.start()
        return self

    def __exit__(self, kind, value, traceback):
        if kind is not None:
            self.cancelled = issubclass(kind, KeyboardInterrupt)
            self.progress.update(self.current, description="Interrupted" if self.cancelled else f"Failed: {value}", visible=True)
        elif self.progress.tasks[self.current].total is None:
            self.progress.update(self.current, description="Finished", total=1, completed=1)
        if kind is None and self.progress.tasks[self.overall].total is None:
            self.progress.update(self.overall, total=1, completed=1)
        self.progress.stop()

    def __call__(self, event):
        key = (event.get("dataset_name", ""), event.get("relative_fov", ""))
        name = "/".join(part for part in key if part)
        kind = event.get("event")
        if kind == "stage":
            task = self.overall if event.get("overall") else self.current
            if self.stage_names.get(task) != event["description"]:
                self.progress.reset(task, total=event.get("total"))
                self.stage_names[task] = event["description"]
            description = event["description"]
            if event.get("overall"):
                description += f" | resumed FOVs {self.resumed} | failed FOVs {self.failed}"
            self.progress.update(task, description=description,
                                 total=event.get("total"), completed=event.get("completed", 0), visible=True)
            return
        if kind == "reconstruction_done":
            if key != self.current_key:
                self.progress.reset(self.current, total=event["total_in_fov"])
                self.current_key = key
            previous = self.counts.get(key, 0)
            completed = event.get("completed_in_fov", previous + 1)
            self.counts[key] = completed
            self.progress.advance(self.overall, max(0, completed - previous))
            self.progress.update(self.current, description=name, total=event["total_in_fov"],
                                 completed=completed, visible=True)
        elif kind == "validation_started":
            self.progress.reset(self.current, total=None)
            self.progress.update(self.current, description=f"{name}: validating outputs", total=None, completed=0, visible=True)
        elif event.get("status") in {"resumed_manifest", "skipped_existing"}:
            count = event.get("reconstruction_count", event.get("movie_count", 0))
            self.progress.advance(self.overall, max(0, count - self.counts.get(key, 0)))
            self.counts[key] = count
            self.resumed += 1
            self.progress.update(self.current, description=f"{name}: resumed", total=count, completed=count, visible=True)
        elif event.get("status") == "failed":
            self.failed += 1
            count = event.get("reconstruction_count", event.get("movie_count", 0))
            self.failed_work += max(0, count - self.counts.get(key, 0))
            stage = event.get("failure_stage", "processing")
            self.progress.update(self.current, description=f"{name}: FAILED ({stage})", visible=True)
            self.progress.console.print(f"FAILED {name} ({stage}): {event.get('error')}", markup=False)
        elif event.get("status") == "written":
            self.progress.update(self.current, description=f"{name}: outputs validated", total=1, completed=1, visible=True)
        self.progress.update(self.overall, description=f"{self.description} | resumed FOVs {self.resumed} | failed FOVs {self.failed} | unfinished jobs {self.failed_work}")


_callback = ContextVar("sacd_progress_callback", default=None)


def report_progress(function):
    """Add an optional callback without passing rendering objects into scientific code."""
    @wraps(function)
    def wrapped(*args, progress_callback=None, **kwargs):
        outermost = _callback.get() is None
        token = _callback.set(progress_callback or _callback.get())
        try:
            return function(*args, **kwargs)
        except BaseException as exc:
            callback = _callback.get()
            if callback is not None and outermost:
                callback({"status": "failed", "relative_fov": function.__name__, "error": str(exc)})
            raise
        finally:
            _callback.reset(token)
    return wrapped


def stage(description, *, completed=0, total=None, overall=False):
    callback = _callback.get()
    if callback is not None:
        callback({"event": "stage", "description": description, "completed": completed,
                  "total": total, "overall": overall})


def report_event(event):
    callback = _callback.get()
    if callback is not None:
        callback(event)


def progress_items(items, description, *, total=None, overall=False):
    if total is None:
        total = len(items)
    stage(description, total=total, overall=overall)
    for index, item in enumerate(items, 1):
        yield item
        stage(description, completed=index, total=total, overall=overall)


def progress_message(*args, **kwargs):
    """Keep CLI messages; notebooks use one updating stage description."""
    update_stage = kwargs.pop("update_stage", True)
    if _callback.get() is None:
        print(*args, **kwargs)
    elif update_stage:
        stage(" ".join(str(arg) for arg in args))
