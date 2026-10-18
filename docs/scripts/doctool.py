# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Helpers for the Claude Project Discipline docs layout.

    uv run docs/scripts/doctool.py lint [--people N] [--strict]
    uv run docs/scripts/doctool.py archive <file.md>... [--dry-run] [--force]
    uv run docs/scripts/doctool.py new decision|spike|plan|handoff|area "<title>"

Layout (relative to the repo root; override the docs dir with --docs):

    docs/architecture/<area>.md              architecture reference
    docs/architecture/decisions/NNNN-*.md    decision records
    docs/architecture/decisions/archive/     superseded / withdrawn decisions
    docs/architecture/decisions/spikes/      spikes
    docs/plans/<topic>.md                    plans
    docs/plans/archive/                      finished plans
    docs/handoffs/<topic>-handoff.md         handoffs

See https://robtaylor.github.io/claude-project-discipline/
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from urllib.parse import quote, unquote

TEMPLATE_NAMES = frozenset({"_template.md", "0000-template.md"})
NON_DOC_NAMES = TEMPLATE_NAMES | {"README.md", "SUMMARY.md"}
ROOT_DOCS = ("CLAUDE.md", "AGENTS.md", "README.md", "CONTRIBUTING.md", "CHANGELOG.md")
DECISION_RE = re.compile(r"^(\d{4})-.+\.md$")

# [text](target "title") and ![alt](target); target may be <angle-bracketed>.
INLINE_LINK_RE = re.compile(r"!?\[(?:[^\[\]]|\[[^\]]*\])*\]\(\s*(<[^>]*>|[^)\s]+)(?:\s+(?:\"[^\"]*\"|'[^']*'))?\s*\)")
REF_DEF_RE = re.compile(r"^[ ]{0,3}\[[^\]]+\]:[ \t]*(<[^>]*>|\S+)", re.MULTILINE)
FENCE_RE = re.compile(r"^[ ]{0,3}(`{3,}|~{3,})")
CODE_SPAN_RE = re.compile(r"(`+)(.+?)\1")
PATH_TOKEN_RE = re.compile(r"(?<![\w/.\-])((?:\.{1,2}/)*[\w.\-]+(?:/[\w.\-]+)*\.md)(?![\w/\-]|\.\w)")
# **Status:** x, **Status (2026-07-06):** x, optionally blockquoted; or a "## Status" section.
STATUS_RE = re.compile(r"^(?:>\s*)?\*\*Status(?:\s*\([^)]*\))?:?\*\*:?\s*(.*)$", re.MULTILINE | re.IGNORECASE)
STATUS_HEADING_RE = re.compile(r"^#{2,}\s*Status\s*$\n+(?:>\s*)?(\S.*)$", re.MULTILINE | re.IGNORECASE)
FIELD_LINE_RE = re.compile(r"^\*\*[^*]+:\*\*")
TABLE_SEP_RE = re.compile(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)*\|?\s*$")


class DocToolError(Exception):
    """A user-facing failure with a message that says what to do next."""


@dataclass(frozen=True)
class DocKind:
    """One kind of doc: where it lives, its template, and its status lifecycle.

    ``statuses`` are the first words a Status line may start with; ``archive_on``
    are the ones that mean "finished, move it to archive/". Keep both in step
    with the template comments and the README status legends.
    """

    name: str
    dir_attr: str
    template: str
    statuses: tuple[str, ...] = ()
    archive_on: tuple[str, ...] = ()
    filename: str = "{slug}.md"


KINDS = {k.name: k for k in (
    DocKind("decision", "decisions", "0000-template.md",
            ("proposed", "accepted", "superseded", "withdrawn"), ("superseded", "withdrawn"), "{num}-{slug}.md"),
    DocKind("spike", "spikes", "_template.md", ("open", "resolved")),
    DocKind("plan", "plans", "_template.md",
            ("proposed", "active", "deferred", "implemented", "withdrawn"), ("implemented", "withdrawn")),
    DocKind("handoff", "handoffs", "_template.md", filename="{slug}-handoff.md"),
    DocKind("area", "architecture", "_template.md"),
)}


@dataclass
class Layout:
    root: Path
    docs: Path
    _texts: dict[Path, str] = field(default_factory=dict, repr=False)

    @property
    def architecture(self) -> Path:
        return self.docs / "architecture"

    @property
    def decisions(self) -> Path:
        return self.architecture / "decisions"

    @property
    def decision_index(self) -> Path:
        return self.decisions / "README.md"

    @property
    def spikes(self) -> Path:
        return self.decisions / "spikes"

    @property
    def plans(self) -> Path:
        return self.docs / "plans"

    @property
    def handoffs(self) -> Path:
        return self.docs / "handoffs"

    def dir_of(self, kind: DocKind) -> Path:
        return getattr(self, kind.dir_attr)

    def markdown_files(self) -> list[Path]:
        """Every .md file under docs/, plus the conventional root-level docs."""
        files = sorted(p for p in self.docs.rglob("*.md") if not _ignored(p, self.docs))
        files += [self.root / n for n in ROOT_DOCS if (self.root / n).is_file()]
        return files

    def text(self, p: Path) -> str:
        """File contents, read once per run."""
        if p not in self._texts:
            self._texts[p] = p.read_text(encoding="utf-8")
        return self._texts[p]

    def rel(self, p: Path) -> str:
        return posix_rel(p, self.root)


def _ignored(p: Path, base: Path) -> bool:
    parts = p.relative_to(base).parts
    return any(part.startswith(".") or part in {"node_modules", "book"} for part in parts)


def norm(p: Path | str) -> Path:
    return Path(os.path.normpath(p))


def posix_rel(target: Path, base: Path) -> str:
    return Path(os.path.relpath(target, base)).as_posix()


def live_docs(d: Path) -> list[Path]:
    """The real docs directly in ``d``: no templates, READMEs, or nav files."""
    return [p for p in sorted(d.glob("*.md")) if p.name not in NON_DOC_NAMES] if d.is_dir() else []


def decision_num(p: Path) -> str:
    m = DECISION_RE.match(p.name)
    assert m, f"{p.name} is not a decision filename"
    return m.group(1)


def decision_files(layout: Layout) -> list[Path]:
    return [p for d in (layout.decisions, layout.decisions / "archive")
            for p in live_docs(d) if DECISION_RE.match(p.name)]


def area_docs(layout: Layout) -> list[Path]:
    return live_docs(layout.architecture)


# --------------------------------------------------------------------------
# Markdown scanning


@dataclass
class Span:
    start: int
    end: int
    text: str


def line_of(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def is_table_row(line: str) -> bool:
    return line.lstrip().startswith("|")


def code_mask(text: str) -> tuple[str, list[bool]]:
    """Blank out fenced blocks and inline code (same length, offsets preserved).

    Also returns, per line, whether that line is inside a fenced block.
    """
    out: list[str] = []
    in_fence: list[bool] = []
    fence: str | None = None
    for line in text.splitlines(keepends=True):
        m = FENCE_RE.match(line)
        if fence is None and not m:
            in_fence.append(False)
            out.append(CODE_SPAN_RE.sub(lambda c: " " * len(c.group(0)), line))
            continue
        in_fence.append(True)
        out.append(re.sub(r"[^\n]", " ", line))
        if fence is None:
            assert m is not None
            fence = m.group(1)
        elif m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence):
            fence = None
    masked = "".join(out)
    assert len(masked) == len(text), "masking must preserve offsets"
    return masked, in_fence


def iter_link_targets(text: str, masked: str | None = None) -> list[Span]:
    """Link targets (inline and reference definitions) outside code."""
    if masked is None:
        masked, _ = code_mask(text)
    spans = []
    for rx in (INLINE_LINK_RE, REF_DEF_RE):
        for m in rx.finditer(masked):
            start, end = m.span(1)
            raw = text[start:end]
            if raw.startswith("<") and raw.endswith(">"):
                start, end, raw = start + 1, end - 1, raw[1:-1]
            spans.append(Span(start, end, raw))
    return sorted(spans, key=lambda s: s.start)


def iter_path_mentions(text: str, skip: list[Span], in_fence: list[bool]) -> list[Span]:
    """`*.md` path tokens in prose and inline code, excluding fenced blocks and link targets."""
    spans = []
    offset = 0
    for line, fenced in zip(text.splitlines(keepends=True), in_fence):
        if not fenced:
            for m in PATH_TOKEN_RE.finditer(line):
                s, e = offset + m.start(1), offset + m.end(1)
                if not any(k.start <= s < k.end for k in skip):
                    spans.append(Span(s, e, m.group(1)))
        offset += len(line)
    return spans


def split_target(target: str) -> tuple[str, str]:
    """('path', '#frag?query') for a link target."""
    m = re.match(r"^([^#?]*)(.*)$", target)
    assert m is not None
    return m.group(1), m.group(2)


def local_target(target: str, base_dir: Path) -> Path | None:
    """Where a link target points, or None if it's external, an anchor, or a `<placeholder>`."""
    if not target or target.startswith(("#", "/")) or "<" in target or ">" in target:
        return None
    if re.match(r"^[a-zA-Z][a-zA-Z0-9+.\-]*:", target):
        return None
    path, _ = split_target(target)
    return norm(base_dir / unquote(path)) if path else None


def local_links(text: str, base_dir: Path) -> Iterator[tuple[Span, Path]]:
    for span in iter_link_targets(text):
        if (target := local_target(span.text, base_dir)) is not None:
            yield span, target


def apply_edits(text: str, edits: list[tuple[int, int, str]]) -> str:
    for start, end, new in sorted(edits, reverse=True):
        text = text[:start] + new + text[end:]
    return text


def status_of(text: str) -> tuple[str | None, int | None]:
    m = STATUS_RE.search(text) or STATUS_HEADING_RE.search(text)
    return (m.group(1).strip(), line_of(text, m.start(1))) if m else (None, None)


def status_word(status: str | None) -> str:
    m = re.match(r"[\W_]*([A-Za-z]+)", status or "")
    return m.group(1).lower() if m else ""


# --------------------------------------------------------------------------
# lint


@dataclass
class Finding:
    path: Path
    line: int
    code: str
    message: str
    severity: str = "error"


@dataclass
class Report:
    layout: Layout
    findings: list[Finding] = field(default_factory=list)

    def add(self, path: Path, line: int | None, code: str, msg: str, severity: str = "error") -> None:
        self.findings.append(Finding(path, line or 1, code, msg, severity))

    def print(self) -> None:
        for f in sorted(self.findings, key=lambda f: (self.layout.rel(f.path), f.line, f.code)):
            print(f"{self.layout.rel(f.path)}:{f.line}: {f.severity} {f.code} {f.message}")


def lint_markdown(layout: Layout, report: Report) -> None:
    """Broken relative links, and table headers with no rows (e.g. a duplicated header)."""
    for path in layout.markdown_files():
        text = layout.text(path)
        for span, target in local_links(text, path.parent):
            if not target.exists():
                report.add(path, line_of(text, span.start), "L001", f"broken link: {span.text}")
        lines = text.splitlines()
        _, in_fence = code_mask(text)
        for i, line in enumerate(lines[:-1]):
            if in_fence[i] or not is_table_row(line) or not TABLE_SEP_RE.match(lines[i + 1]):
                continue
            if not (i + 2 < len(lines) and is_table_row(lines[i + 2])):
                report.add(path, i + 1, "M001", "table header has no rows (duplicated or orphaned header?)", "warning")


def lint_statuses(layout: Layout, report: Report) -> None:
    """Every doc of a kind with a lifecycle has a known Status, and sits in archive/ iff it's finished."""
    for kind in KINDS.values():
        if not kind.statuses:
            continue
        home = layout.dir_of(kind)
        for archived, d in ((False, home), (True, home / "archive")):
            for path in live_docs(d):
                if kind.name == "decision" and not DECISION_RE.match(path.name):
                    continue
                status, sline = status_of(layout.text(path))
                word = status_word(status)
                if status is None:
                    report.add(path, 1, "S001", f"{kind.name} has no **Status:** line")
                    continue
                if word not in kind.statuses:
                    report.add(path, sline, "S001", f"{kind.name} status '{status}' doesn't start with one of: "
                               f"{', '.join(s.title() for s in kind.statuses)}", "warning")
                if not kind.archive_on:
                    continue
                if word in kind.archive_on and not archived:
                    report.add(path, sline, "R001",
                               f"status is {word.title()}; move it with: doctool.py archive {layout.rel(path)}")
                elif archived and word not in kind.archive_on:
                    report.add(path, sline, "R001", f"in archive/ but status is not "
                               f"{' or '.join(s.title() for s in kind.archive_on)}", "warning")


def lint_decisions(layout: Layout, report: Report) -> None:
    index = layout.decision_index
    index_rows = parse_index(layout.text(index), layout.decisions) if index.is_file() else {}
    areas = area_docs(layout)
    linked_from_areas = {t for area in areas for _, t in local_links(layout.text(area), area.parent)}
    seen: dict[str, Path] = {}

    for path in decision_files(layout):
        text = layout.text(path)
        num = decision_num(path)
        if num in seen:
            report.add(path, 1, "D001", f"decision number {num} also used by {layout.rel(seen[num])}")
        seen[num] = path
        status, _ = status_of(text)
        word = status_word(status)

        for m in re.finditer(r"^(?:#+\s*Amendment\b|>\s*\*\*Amendment\b)", text, re.MULTILINE | re.IGNORECASE):
            report.add(path, line_of(text, m.start()), "D004", "dated Amendment block: decisions are immutable; "
                       "record the change in the area doc or a new decision", "warning")

        if path.parent.name != "archive":
            lint_current_architecture(report, path, text, word)
            if word == "accepted" and areas and path not in linked_from_areas:
                report.add(path, 1, "A002", "accepted decision is not linked from any area doc's "
                           "'Decisions behind this'", "warning")

        row = index_rows.get(path)
        if index.is_file() and row is None:
            report.add(index, 1, "D002", f"decision {num} ({layout.rel(path)}) is missing from the index")
        elif row is not None and status_word(row[1]) and word and status_word(row[1]) != word:
            report.add(index, row[0], "D002", f"index says '{row[1]}' but {path.name} says '{status}'")


def lint_current_architecture(report: Report, path: Path, text: str, word: str) -> None:
    lines = text.splitlines()
    idx = next((i for i, ln in enumerate(lines) if ln.startswith("**Current architecture:**")), None)
    if idx is None:
        severity = "warning" if word == "proposed" else "error"
        report.add(path, 1, "D005", "no **Current architecture:** forward-link to an area doc", severity)
        return
    # The forward-link is its own one-line paragraph. If prose follows it with no
    # blank line, it was inserted into the middle of another paragraph.
    if idx + 1 < len(lines) and lines[idx + 1].strip() and not FIELD_LINE_RE.match(lines[idx + 1]):
        report.add(path, idx + 2, "D006", "prose continues after the **Current architecture:** line; "
                   "was it inserted mid-paragraph (e.g. splitting the Status)?")


def parse_index(text: str, base: Path) -> dict[Path, tuple[int, str]]:
    """Map each decision linked from an index table row to (line, status cell)."""
    rows: dict[Path, tuple[int, str]] = {}
    for lineno, line in enumerate(text.splitlines(), 1):
        if not is_table_row(line) or TABLE_SEP_RE.match(line):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        for _, target in local_links(line, base):
            if DECISION_RE.match(target.name):
                rows.setdefault(target, (lineno, cells[-1] if len(cells) >= 3 else ""))
    return rows


def lint_handoffs(layout: Layout, report: Report, people: int | None) -> None:
    live = live_docs(layout.handoffs)
    for p in live:
        if not p.name.endswith("-handoff.md"):
            report.add(p, 1, "H002", "handoff file should be named <topic>-handoff.md", "warning")
    if people is not None and len(live) > 2 * people:
        report.add(layout.handoffs, 1, "H001",
                   f"{len(live)} live handoffs for {people} people (> 2x): resolve, consolidate, or drop threads",
                   "warning")


def lint_areas(layout: Layout, report: Report) -> None:
    for path in area_docs(layout):
        text = layout.text(path)
        for heading in ("Implementation status", "Decisions behind this"):
            if not re.search(rf"^##\s+{heading}\b", text, re.MULTILINE | re.IGNORECASE):
                report.add(path, 1, "A001", f"area doc has no '## {heading}' section", "warning")


def cmd_lint(layout: Layout, args: argparse.Namespace) -> int:
    report = Report(layout)
    lint_markdown(layout, report)
    lint_statuses(layout, report)
    lint_decisions(layout, report)
    lint_handoffs(layout, report, args.people)
    lint_areas(layout, report)
    report.print()
    errors = sum(f.severity == "error" for f in report.findings)
    warnings = len(report.findings) - errors
    print(f"{errors} error(s), {warnings} warning(s)", file=sys.stderr)
    return 1 if errors or (args.strict and warnings) else 0


# --------------------------------------------------------------------------
# archive


def tracked_text_files_naming(layout: Layout, names: list[str]) -> set[Path]:
    """Tracked non-binary files that mention any of ``names`` (cheap pre-filter)."""
    args = [a for n in names for a in ("-e", n)]
    r = subprocess.run(["git", "grep", "-lzIF", *args], check=False, cwd=layout.root,
                       capture_output=True, text=True)
    return {norm(layout.root / p) for p in r.stdout.split("\0") if p}


def root_path_re(rel: str) -> re.Pattern[str]:
    """A repo-root path like docs/plans/x.md, not as part of a longer path."""
    return re.compile(rf"(?:(?<![\w\-./])|(?<=\./)){re.escape(rel)}(?![\w\-/]|\.\w)")


def rewrite_md(layout: Layout, text: str, old_dir: Path, new_dir: Path, moves: dict[Path, Path]) -> tuple[str, int]:
    """Repoint links and path mentions in a Markdown file that lived in ``old_dir``
    and now lives in ``new_dir``. Every path is resolved against ``old_dir``,
    before anything moves, so a same-named file at another depth can't be hit."""
    masked, in_fence = code_mask(text)
    links = iter_link_targets(text, masked)
    edits = []
    for span in links:
        target = local_target(span.text, old_dir)
        if target is None:
            continue
        new_target = moves.get(target, target)
        if new_target == target and new_dir == old_dir:
            continue
        new = posix_rel(new_target, new_dir) + split_target(span.text)[1]
        if new != span.text:
            edits.append((span.start, span.end, new))
    for span in iter_path_mentions(text, links, in_fence):
        rel_target = norm(old_dir / span.text)
        root_target = norm(layout.root / span.text)
        if root_target in moves:
            continue  # repo-root style; handled by the root-path pass
        if rel_target in moves:
            edits.append((span.start, span.end, posix_rel(moves[rel_target], new_dir)))
        elif new_dir != old_dir and rel_target.exists() and rel_target != root_target:
            edits.append((span.start, span.end, posix_rel(rel_target, new_dir)))
    return apply_edits(text, edits), len(edits)


def plan_archive(layout: Layout, moves: dict[Path, Path]) -> tuple[dict[Path, str], list[str], list[str]]:
    """Work out every file edit for ``moves`` (src -> dest) without touching disk.

    Returns (new text keyed by final path, notes, leftover mentions to check by hand).
    """
    names = sorted({n for s in moves for n in (s.name, quote(s.name))})
    candidates = set(layout.markdown_files()) | set(moves) | tracked_text_files_naming(layout, names)
    root_res = [(root_path_re(layout.rel(s)), layout.rel(d)) for s, d in moves.items()]
    leftover_re = re.compile("|".join(rf"\b{re.escape(s.parent.name)}/{re.escape(s.name)}(?![\w\-/]|\.\w)"
                                      for s in moves))

    changes: dict[Path, str] = {}
    notes: list[str] = []
    leftovers: list[str] = []
    for path in sorted(candidates):
        if not path.is_file():
            continue
        text = layout.text(path)
        if path not in moves and not any(n in text for n in names):
            continue  # can't refer to a moved file without naming it
        final = moves.get(path, path)
        new, count = text, 0
        if path.suffix == ".md":
            new, count = rewrite_md(layout, text, path.parent, final.parent, moves)
            stale = [t for _, t in local_links(new, final.parent) if t in moves]
            assert not stale, f"{layout.rel(final)} still links to an old path: {stale}"
        for rx, new_rel in root_res:
            new, n = rx.subn(new_rel, new)
            count += n
        if new != text:
            changes[final] = new
            what = "outbound reference(s) re-depthed" if path in moves else "reference(s) repointed"
            notes.append(f"  {layout.rel(final)}: {count} {what}")
        for lineno, line in enumerate(new.splitlines(), 1):
            if leftover_re.search(line):
                leftovers.append(f"  {layout.rel(final)}:{lineno}: {line.strip()}")
    return changes, notes, leftovers


def check_archivable(layout: Layout, src: Path, force: bool) -> None:
    kind = next((k for k in KINDS.values() if k.archive_on and src.parent == norm(layout.dir_of(k))), None)
    if kind is None or (kind.name == "decision" and not DECISION_RE.match(src.name)):
        homes = " or ".join(f"{layout.rel(layout.dir_of(k))}/" for k in KINDS.values() if k.archive_on)
        raise DocToolError(f"{layout.rel(src)} is not a doc that can be archived (expected one in {homes})")
    status, _ = status_of(layout.text(src))
    if status_word(status) not in kind.archive_on and not force:
        raise DocToolError(
            f"{layout.rel(src)} has status '{status}'. Set it to one of "
            f"{', '.join(a.title() for a in kind.archive_on)} first, or pass --force.")


def cmd_archive(layout: Layout, args: argparse.Namespace) -> int:
    moves: dict[Path, Path] = {}
    for arg in args.files:
        src = norm(Path(arg).absolute())
        dest = src.parent / "archive" / src.name
        if src.parent.name == "archive" or (not src.exists() and dest.exists()):
            print(f"skip: {arg} is already archived")
            continue
        if not src.is_file():
            raise DocToolError(f"{arg}: no such file")
        if dest.exists():
            raise DocToolError(f"{layout.rel(dest)} already exists")
        check_archivable(layout, src, args.force)
        moves[src] = dest
    if not moves:
        return 0

    changes, notes, leftovers = plan_archive(layout, moves)
    for src, dest in moves.items():
        print(f"{'would move' if args.dry_run else 'move'}: {layout.rel(src)} -> {layout.rel(dest)}")
    print("\n".join(notes) if notes else "  no references to repoint")
    if not args.dry_run:
        for src, dest in moves.items():
            dest.parent.mkdir(parents=True, exist_ok=True)
            git_mv(layout, src, dest)
        for path, text in changes.items():
            path.write_text(text, encoding="utf-8")
    if leftovers:
        print(f"\n{len(leftovers)} mention(s) left for a manual check:")
        print("\n".join(leftovers))
    if not args.dry_run:
        print("\nDone. Regroup index tables / SUMMARY.md headings for finished items, then run lint.")
    return 0


def git_mv(layout: Layout, src: Path, dest: Path) -> None:
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(src)], check=False,
                             cwd=layout.root, capture_output=True).returncode == 0
    if tracked:
        subprocess.run(["git", "mv", str(src), str(dest)], cwd=layout.root, check=True)
    else:
        src.rename(dest)


# --------------------------------------------------------------------------
# new


def slugify(title: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
    if not slug:
        raise DocToolError(f"can't make a filename from title {title!r}")
    return slug


def cmd_new(layout: Layout, args: argparse.Namespace) -> int:
    kind = KINDS[args.kind]
    home = layout.dir_of(kind)
    template = home / kind.template
    is_decision = kind.name == "decision"
    num = f"{max((int(decision_num(p)) for p in decision_files(layout)), default=0) + 1:04d}" if is_decision else ""
    dest = home / kind.filename.format(slug=slugify(args.title), num=num)
    if not template.is_file():
        raise DocToolError(f"template {layout.rel(template)} not found")
    if dest.exists():
        raise DocToolError(f"{layout.rel(dest)} already exists")

    lines = layout.text(template).splitlines(keepends=True)
    assert lines and lines[0].startswith("# "), f"{template} should start with an H1"
    if is_decision:
        lines[0] = f"# Decision {num} — {args.title}\n"
    else:
        lines[0] = re.sub(r"<[^>]+>", args.title, lines[0], count=1)
    text = "".join(lines)
    if kind.name == "handoff":
        text = text.replace("YYYY-MM-DD", date.today().isoformat(), 1)
    dest.write_text(text, encoding="utf-8")
    print(f"Created {layout.rel(dest)}")
    if is_decision:
        add_index_row(layout, num, dest, args.title)
    return 0


def add_index_row(layout: Layout, num: str, dest: Path, title: str) -> None:
    index = layout.decision_index
    if not index.is_file():
        print(f"No {layout.rel(index)}; add the index row by hand.")
        return
    lines = layout.text(index).splitlines(keepends=True)
    start = next((i for i, ln in enumerate(lines) if re.match(r"^##\s+Index\b", ln)), None)
    last = None
    for i in range(start + 1 if start is not None else len(lines), len(lines)):
        if is_table_row(lines[i]):
            last = i
        elif (last is not None and lines[i].strip()) or lines[i].startswith("#"):
            break
    if last is None:
        print(f"No table under '## Index' in {layout.rel(index)}; add the row by hand.")
        return
    lines.insert(last + 1, f"| [{num}]({dest.name}) | {title} | Proposed |\n")
    index.write_text("".join(lines), encoding="utf-8")
    print(f"Added index row to {layout.rel(index)}")


# --------------------------------------------------------------------------


def find_root(start: Path) -> Path:
    r = subprocess.run(["git", "rev-parse", "--show-toplevel"], check=False, cwd=start,
                       capture_output=True, text=True)
    return Path(r.stdout.strip()) if r.returncode == 0 else start


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, help="repo root (default: git toplevel, else cwd)")
    ap.add_argument("--docs", default="docs", help="docs directory relative to the root (default: docs)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    lp = sub.add_parser("lint", help="check the docs against the discipline's rules")
    lp.add_argument("--people", type=int, help="people actively working, for the handoff WIP check")
    lp.add_argument("--strict", action="store_true", help="fail on warnings too")

    arp = sub.add_parser("archive", help="move finished plans or superseded decisions to archive/ and repoint links")
    arp.add_argument("files", nargs="+", metavar="FILE")
    arp.add_argument("--dry-run", action="store_true")
    arp.add_argument("--force", action="store_true", help="archive even if the status doesn't say it's finished")

    np_ = sub.add_parser("new", help="create a doc from its template")
    np_.add_argument("kind", choices=list(KINDS))
    np_.add_argument("title")

    args = ap.parse_args(argv)
    root = norm((args.root or find_root(Path.cwd())).absolute())
    layout = Layout(root=root, docs=norm(root / args.docs))
    if not layout.docs.is_dir():
        print(f"error: {layout.docs} is not a directory", file=sys.stderr)
        return 2
    try:
        return {"lint": cmd_lint, "archive": cmd_archive, "new": cmd_new}[args.cmd](layout, args)
    except DocToolError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
