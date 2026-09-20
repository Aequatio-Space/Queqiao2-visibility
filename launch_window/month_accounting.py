"""launch_window/month_accounting.py — C-8 calendar-month accounting (PRD v6 §9).

F-8 / F-11 fixes:

- ``n_days_in_month`` is computed from the Gregorian calendar, never hardcoded
  and never derived from ``(end_et - start_et) / 86400``.
- Coverage is ``月内窗口日数 ÷ 月内天数``.  Window days that fall outside the
  calendar month (the neighbouring cluster) are reported separately as
  ``adjacent_cluster_days`` and are **not** counted in the denominator.
- A single source of truth for the calendar-month span strings so every scan
  aligns to ``YYYY-MM-01T00:00:00Z → YYYY-MM+1-01T00:00:00Z``.

This module is deliberately free of SPICE/geometry imports so it can be
unit-tested without loading kernels.
"""
from __future__ import annotations

import calendar
from datetime import datetime, timezone


def days_in_month(year: int, month: int) -> int:
    """Number of days in a Gregorian calendar month (no hardcoded table)."""
    if not (1 <= month <= 12):
        raise ValueError(f"month must be 1..12, got {month}")
    return calendar.monthrange(year, month)[1]


def month_start_utc(year: int, month: int) -> str:
    """ISO UTC string for the first instant of a calendar month."""
    if not (1 <= month <= 12):
        raise ValueError(f"month must be 1..12, got {month}")
    return f"{year:04d}-{month:02d}-01T00:00:00Z"


def next_month_start_utc(year: int, month: int) -> str:
    """ISO UTC string for the first instant of the following calendar month."""
    if not (1 <= month <= 12):
        raise ValueError(f"month must be 1..12, got {month}")
    if month == 12:
        return f"{year + 1:04d}-01-01T00:00:00Z"
    return f"{year:04d}-{month + 1:02d}-01T00:00:00Z"


def month_prefix(year: int, month: int) -> str:
    """Date prefix ``YYYY-MM-`` used to classify window-day strings."""
    return f"{year:04d}-{month:02d}-"


def split_window_days(
    window_days: list[str], year: int, month: int
) -> tuple[list[str], list[str]]:
    """Split window-day strings into (in_month, adjacent) lists.

    Adjacent days are those belonging to the neighbouring cluster that the
    scan straddles; per C-8 they must not enter the month's denominator.
    """
    prefix = month_prefix(year, month)
    in_month = [d for d in window_days if d.startswith(prefix)]
    adjacent = [d for d in window_days if not d.startswith(prefix)]
    return sorted(in_month), sorted(adjacent)


def monthly_coverage(
    window_days: list[str], year: int, month: int
) -> dict:
    """C-8 coverage accounting for a single calendar month.

    Returns ``n_days_in_month`` (calendar-derived), the in-month window days,
    the adjacent cluster days, and ``coverage = in_month / n_days_in_month``.
    """
    n_days = days_in_month(year, month)
    in_month, adjacent = split_window_days(window_days, year, month)
    coverage = round(len(in_month) / n_days, 4) if n_days else None
    return {
        "year": year,
        "month": month,
        "n_days_in_month": n_days,
        "n_window_days_in_month": len(in_month),
        "window_days_in_month": in_month,
        "adjacent_cluster_days": adjacent,
        "coverage": coverage,
    }


def coverage_recount_diff(
    per_mission: list[dict], before_key: str, after_key: str
) -> dict:
    """C-8 before/after coverage recount.

    ``per_mission`` entries carry the frozen value under ``before_key`` and
    the freshly recomputed value under ``after_key``.  Returns a map
    ``mission -> {before, after, delta}`` plus a summary ``n_changed``.
    ``None`` before-values (e.g. a mission absent from the frozen baseline)
    are reported with a ``None`` delta rather than a fabricated zero.
    """
    diff = {}
    n_changed = 0
    for m in per_mission:
        before_raw = m.get(before_key)
        after_raw = m.get(after_key)
        before = round(float(before_raw), 4) if before_raw is not None else None
        after = round(float(after_raw), 4) if after_raw is not None else None
        if before is None and after is None:
            delta = None
        else:
            delta = round(after - before, 4) if (before is not None and after is not None) else None
        if delta is not None and delta != 0.0:
            n_changed += 1
        diff[m["mission"]] = {"before": before, "after": after, "delta": delta}
    return {"per_mission": diff, "n_changed": n_changed}


def now_utc_iso() -> str:
    """Current UTC instant as an ISO-8601 second-precision string."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
