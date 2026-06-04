import subprocess
from pathlib import Path


def get_git_commit():
    """
    Return current git commit hash.
    Works for editable install.
    """

    try:
        repo_path = Path(__file__).resolve().parents[2]

        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_path,
        ).decode().strip()

        return commit

    except Exception:
        return "unknown"