"""Ticket-driven orchestrator for the agents.md workflow used in bwa-rs.

A ticket is considered done only when its PR (titled "T{N}: ...") is merged
on the remote, with CI checks passing. The agent runs as a fresh subprocess
per ticket so its context is destroyed on exit.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

STATUS_PENDING = "⬜ Pending"
STATUS_IN_PROGRESS = "🟡 In Progress"
STATUS_DONE = "✅ Done"
STATUS_FAILED = "❌ Failed"

TICKET_RE = re.compile(
    r"^### T(?P<n>\d+):\s*(?P<title>.+?)\s*\n"
    r"(?:\s*\n)*"
    r"\*\*Status:\*\*\s*(?P<status>.+?)\s*\n",
    re.MULTILINE,
)


@dataclass
class Ticket:
    """A ticket header parsed from tickets.md.

    Attributes:
        number: Integer ticket number (e.g. 29 for T29).
        title: Ticket title text.
        status: Current status string including emoji.
    """

    number: int
    title: str
    status: str

    @property
    def tag(self) -> str:
        """Returns the canonical T{N} tag."""
        return f"T{self.number}"


def parse_tickets(path: Path) -> list[Ticket]:
    """Parses ticket headers from tickets.md.

    Args:
        path: Path to tickets.md.

    Returns:
        Tickets in document order.
    """
    text = path.read_text()
    return [
        Ticket(int(m.group("n")), m.group("title").strip(), m.group("status").strip())
        for m in TICKET_RE.finditer(text)
    ]


def next_pending(tickets: list[Ticket]) -> Ticket | None:
    """Returns the first pending ticket, or None if the queue is empty."""
    for t in tickets:
        if "Pending" in t.status:
            return t
    return None


def update_status(path: Path, ticket: Ticket, new_status: str) -> None:
    """Updates a single ticket's Status line in place.

    Args:
        path: Path to tickets.md.
        ticket: Ticket to update.
        new_status: Replacement status string (with emoji).
    """
    text = path.read_text()
    pattern = re.compile(
        rf"(### T{ticket.number}:[^\n]*\n(?:[ \t]*\n)*\*\*Status:\*\*[ \t]*)[^\n]+",
    )
    new_text, n = pattern.subn(rf"\g<1>{new_status}", text, count=1)
    if n != 1:
        raise ValueError(f"Could not update status for T{ticket.number}")
    path.write_text(new_text)


def run_agent(ticket: Ticket, agent_cmd: list[str], workdir: Path, timeout: int) -> int:
    """Spawns the coding agent for a single ticket and waits for exit.

    The prompt directs the agent to follow agents.md, which prescribes the full
    branch/test/implement/simplify/PR loop.

    Args:
        ticket: The ticket being worked on.
        agent_cmd: Command list to invoke the coding agent.
        workdir: Repository root.
        timeout: Maximum wall-clock seconds.

    Returns:
        Subprocess exit code; -1 on timeout.
    """
    reset = subprocess.run(
        ["git", "checkout", "master"], cwd=workdir, capture_output=True, text=True
    )
    if reset.returncode != 0:
        print(f"[runner] git checkout master failed: {reset.stderr}", file=sys.stderr)
        return reset.returncode
    subprocess.run(["git", "pull", "--ff-only"], cwd=workdir, check=False)

    prompt = (
        f"Pull ticket {ticket.tag} ({ticket.title}) following agents.md. "
        "Plan, execute, then simplify before creating the PR. "
        "CI must pass before merge. Only spawn one tool at a time. "
        "CRITICAL: emit exactly one tool call per response, never multiple in parallel."
    )
    print(f"[runner] launching agent for {ticket.tag}", file=sys.stderr)
    try:
        proc = subprocess.run(
            [*agent_cmd, prompt], capture_output=True, text=True, cwd=workdir,
            timeout=timeout, check=False,
        )
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        if proc.returncode != 0 and "Failed to parse input" in (proc.stdout + proc.stderr):
            return -2
        return proc.returncode
    except subprocess.TimeoutExpired:
        print(f"[runner] {ticket.tag} timed out after {timeout}s", file=sys.stderr)
        return -1


def find_pr(ticket: Ticket, workdir: Path) -> dict | None:
    """Looks up the PR whose title starts with this ticket's tag.

    Args:
        ticket: The ticket to search for.
        workdir: Repository root.

    Returns:
        Parsed gh PR JSON (number, state, mergedAt, etc.), or None.
    """
    result = subprocess.run(
        ["gh", "pr", "list", "--state", "all",
         "--search", f"{ticket.tag}: in:title",
         "--json", "number,title,state,mergedAt,headRefName",
         "--limit", "5"],
        cwd=workdir, capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return None
    prs = json.loads(result.stdout or "[]")
    for pr in prs:
        if pr["title"].startswith(f"{ticket.tag}:"):
            return pr
    return None


def wait_for_merge(
    ticket: Ticket, workdir: Path, poll_interval: int, max_wait: int
) -> str:
    """Polls until the ticket's PR is merged, fails, or the deadline passes.

    Args:
        ticket: Ticket whose PR to monitor.
        workdir: Repository root.
        poll_interval: Seconds between polls.
        max_wait: Maximum total seconds to wait.

    Returns:
        One of "merged", "closed", "ci_failed", "no_pr", "timeout".
    """
    deadline = time.monotonic() + max_wait
    last_seen_pr = None
    while time.monotonic() < deadline:
        pr = find_pr(ticket, workdir)
        if pr is None:
            time.sleep(poll_interval)
            continue
        last_seen_pr = pr
        if pr["state"] == "MERGED" or pr.get("mergedAt"):
            return "merged"
        if pr["state"] == "CLOSED":
            return "closed"
        checks = subprocess.run(
            ["gh", "pr", "checks", str(pr["number"]), "--json", "state,conclusion"],
            cwd=workdir, capture_output=True, text=True, check=False,
        )
        if checks.returncode == 0 and checks.stdout.strip():
            try:
                states = json.loads(checks.stdout)
                if any(c.get("conclusion") in {"FAILURE", "CANCELLED"} for c in states):
                    return "ci_failed"
            except json.JSONDecodeError:
                pass
        time.sleep(poll_interval)
    return "no_pr" if last_seen_pr is None else "timeout"


def process_one(
    path: Path, agent_cmd: list[str], workdir: Path,
    agent_timeout: int, ci_timeout: int, poll_interval: int,
) -> dict:
    """Processes a single pending ticket through the full agents.md loop.

    Args:
        path: Path to tickets.md.
        agent_cmd: Coding agent invocation.
        workdir: Repository root.
        agent_timeout: Seconds to allow agent subprocess to run.
        ci_timeout: Seconds to wait for PR merge after agent exits.
        poll_interval: Seconds between gh polls.

    Returns:
        Dict with ticket tag and outcome, or {"done": True} when queue empty.
    """
    tickets = parse_tickets(path)
    ticket = next_pending(tickets)
    if ticket is None:
        return {"done": True}

    update_status(path, ticket, STATUS_IN_PROGRESS)
    rc = -2
    attempts = 0
    while rc == -2 and attempts < 3:
        attempts += 1
        if attempts > 1:
            print(f"[runner] retry {attempts}/3 for {ticket.tag} after parse error", file=sys.stderr)
        rc = run_agent(ticket, agent_cmd, workdir, agent_timeout)
    if rc != 0:
        update_status(path, ticket, STATUS_FAILED)
        reason = "parse_error_max_retries" if rc == -2 else "agent_exit"
        return {"ticket": ticket.tag, "outcome": reason, "rc": rc}

    outcome = wait_for_merge(ticket, workdir, poll_interval, ci_timeout)
    final = STATUS_DONE if outcome == "merged" else STATUS_FAILED
    update_status(path, ticket, final)
    return {"ticket": ticket.tag, "outcome": outcome}


def main(
    path: Path, agent_cmd: list[str], workdir: Path,
    agent_timeout: int, ci_timeout: int, poll_interval: int, halt_on_fail: bool,
) -> None:
    """Drives the queue until empty or a ticket fails.

    Args:
        path: tickets.md path.
        agent_cmd: Coding agent invocation.
        workdir: Repository root.
        agent_timeout: Per-ticket agent timeout in seconds.
        ci_timeout: Per-ticket CI/merge timeout in seconds.
        poll_interval: Polling cadence in seconds.
        halt_on_fail: If True, stop on first non-merged ticket.
    """
    while True:
        result = process_one(
            path, agent_cmd, workdir, agent_timeout, ci_timeout, poll_interval,
        )
        print(json.dumps(result), flush=True)
        if result.get("done"):
            return
        if halt_on_fail and result["outcome"] != "merged":
            return


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--tickets", type=Path, required=True)
    p.add_argument("--workdir", type=Path, default=Path.cwd())
    p.add_argument("--agent-timeout", type=int, default=3600)
    p.add_argument("--ci-timeout", type=int, default=1800)
    p.add_argument("--poll-interval", type=int, default=20)
    p.add_argument("--no-halt", action="store_true")
    p.add_argument("agent_cmd", nargs=argparse.REMAINDER)
    args = p.parse_args()
    main(
        args.tickets, args.agent_cmd, args.workdir,
        args.agent_timeout, args.ci_timeout, args.poll_interval,
        halt_on_fail=not args.no_halt,
    )
