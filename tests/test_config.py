from pathlib import Path
import os

from blr.cli.config import make_paths_absolute
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
