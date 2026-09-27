"""Collect what went wrong across a batch, so a run asks once.

Stopping at the first bad sequence in a batch of thirty wastes the other
twenty-nine. Every skill gathers its problems here instead, prints one table at
the end, and asks once.
"""


class Report:
    def __init__(self):
        self.blocked = []
        self.warned = []

    def block(self, name, reason):
        """Something that has to be settled before this sequence can go on."""
        self.blocked.append((name, reason))

    def warn(self, name, reason):
        """Something worth knowing that does not stop the work."""
        self.warned.append((name, reason))

    def add_warnings(self, name, reasons):
        for reason in reasons:
            self.warn(name, reason)

    def table(self):
        rows = ([(name, "stops", why) for name, why in self.blocked]
                + [(name, "warning", why) for name, why in self.warned])
        if not rows:
            return ""
        width = max(len(name) for name, _, _ in rows)
        lines = [f"{'sequence':{width}s}  {'kind':8s}  what"]
        lines += [f"{name:{width}s}  {kind:8s}  {why}" for name, kind, why in rows]
        return "\n".join(lines)

    def question(self):
        """The one thing to ask, or an empty string when nothing is blocking."""
        if not self.blocked:
            return ""
        names = ", ".join(name for name, _ in self.blocked)
        return (f"{len(self.blocked)} of these need a decision before I go on: "
                f"{names}. How would you like to handle them?")
