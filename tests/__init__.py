"""Empty on purpose, and needed.

Without it, `python3 -m unittest discover -s tests -t .` stops with "Start
directory is not importable" before running anything. That command is the only
way to run all five suites in one go, and it is the one pre-allowlisted in
.claude/settings.local.json, so this file is what makes it work.

Nothing imports from here. Every suite is also runnable on its own, the way
README.md lists them.
"""
