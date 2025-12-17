"""Utility functions for report generation."""


def convert_bp(size):
    """Convert bp values to appropriate metric prefix."""
    size = float(str(size).replace(",", ""))
    for x in ["bp", "Kb", "Mb", "Gb", "Tb"]:
        if size < 1000.0:
            if x == "bp":
                return "{:.0f} {}".format(size, x)
            else:
                return "{:.2f} {}".format(size, x)
        size /= 1000.0
    return size


def fill_none(val, fill_na):
    """Return fill_na instead of None fields for display."""
    return val if val is not None else fill_na


def capitalize_name(name):
    """Capitalize each word in the name."""
    if not name:
        return name
    words = name.split()
    capitalized = []
    for word in words:
        if len(word) == 1:
            capitalized.append(word.upper())
        else:
            capitalized.append(word.capitalize())
    return ' '.join(capitalized)
