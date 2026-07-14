import os
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from meew_annotation_paths import expand_path, resolve_resource_path


class ResourceResolutionTests(unittest.TestCase):
    def test_resource_precedence_and_paths_with_spaces(self):
        with patch.dict(
            os.environ,
            {
                "HOME": "/tmp/annotation home",
                "XDG_DATA_HOME": "/tmp/annotation xdg",
                "MEEW_RESOURCES": "/tmp/annotation resources",
                "GTDBTK_DATA_PATH": "/tmp/tool GTDB",
            },
            clear=True,
        ):
            self.assertEqual(
                resolve_resource_path(
                    "/tmp/explicit GTDB", "GTDBTK_DATA_PATH", "gtdbtk_db"
                ),
                "/tmp/explicit GTDB",
            )
            self.assertEqual(
                resolve_resource_path("", "GTDBTK_DATA_PATH", "gtdbtk_db"),
                "/tmp/tool GTDB",
            )
            del os.environ["GTDBTK_DATA_PATH"]
            self.assertEqual(
                resolve_resource_path("", "GTDBTK_DATA_PATH", "gtdbtk_db"),
                "/tmp/annotation resources/gtdbtk_db",
            )
            del os.environ["MEEW_RESOURCES"]
            self.assertEqual(
                resolve_resource_path("", "GTDBTK_DATA_PATH", "gtdbtk_db"),
                "/tmp/annotation xdg/meew/gtdbtk_db",
            )
            del os.environ["XDG_DATA_HOME"]
            self.assertEqual(
                resolve_resource_path("", "GTDBTK_DATA_PATH", "gtdbtk_db"),
                "/tmp/annotation home/.local/share/meew/gtdbtk_db",
            )
            self.assertEqual(
                expand_path("relative bins/input.fa"), "relative bins/input.fa"
            )


if __name__ == "__main__":
    unittest.main()
