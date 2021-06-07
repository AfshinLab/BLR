from pathlib import Path
import os

from blr.cli.config import make_paths_absolute, flatten, update_changes_set
from .utils import tempinput


def test_make_paths_absolute():
    absolute_path = str(Path(__file__))
    absolute_unresolved_path = str(Path(__file__).parent / ".." / Path(__file__).parents[0] / Path(__file__).name)
    relative_path = str(".." / Path(__file__).parents[0] / Path(__file__).name)
    not_a_path = "test.txt"
    assert make_paths_absolute(absolute_path) == absolute_path
    assert make_paths_absolute(absolute_unresolved_path) == absolute_path
    assert make_paths_absolute(relative_path) == absolute_path
    assert make_paths_absolute(not_a_path) == not_a_path

    # Don't resolve symlinks
    with tempinput(b" ") as file:
        file = Path(file)
        symlink_path = str(file.parent / "tmp_symlink")
        os.symlink(file, symlink_path)
        relative_symlink_path = str(file.parent / ".." / file.parents[0] / "tmp_symlink")
        assert make_paths_absolute(relative_symlink_path) == symlink_path
        os.unlink(symlink_path)


def test_flatten():
    nested = {"key": {"subkey": "value"}}
    flattened = flatten(nested)
    assert "key.subkey" in flattened
    assert flattened["key.subkey"] == "value"


def test_update_changes_set():
    configs = {
        "foo": 3,
        "bar": 4,
        "nested": {"key": "new_value"}
    }
    changes_set = [("bar", "1"), ("nested.key", "old_value")]
    changes_set = update_changes_set(changes_set, configs)

    changes = dict(changes_set)
    assert changes["foo"] == "3"
    assert changes["bar"] == "1"
    assert changes["nested.key"] == "old_value"
