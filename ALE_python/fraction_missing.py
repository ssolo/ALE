"""Utility for reading ALE fraction missing files."""


def read_fraction_missing_file(filepath: str) -> dict[str, float]:
    """Read a fraction missing file and return a dict mapping species to fraction values.

    The file format is one line per species:
        species_name:fraction_value

    Args:
        filepath: Path to the fraction missing file. If empty string, returns empty dict.

    Returns:
        A dict mapping species name (str) to fraction value (float between 0 and 1).
    """
    if filepath == "":
        return {}

    result = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            species_name, fraction_value = line.split(":")
            result[species_name] = float(fraction_value)

    return result
