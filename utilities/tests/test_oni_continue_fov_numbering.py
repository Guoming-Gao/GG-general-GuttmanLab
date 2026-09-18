import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    'oni', Path(__file__).resolve().parents[1] / 'oni_continue_fov_numbering.py')
oni = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(oni)


class NumberingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.before = self.root / 'earlier date'
        self.after = self.root / 'later date'
        self.before.mkdir()
        self.after.mkdir()

    def fov(self, root, name, contents=False):
        path = root / name
        path.mkdir()
        if contents:
            nested = path / ('nested_' + name)
            nested.mkdir()
            (nested / ('experiment_' + name + '_pos0.tif')).write_bytes(b'raw\x00data')
            (nested / 'metadata.json').write_text('{"original": "' + name + '"}')
        return path

    def plan(self):
        return oni.build_plan(self.before, self.after)

    def test_recursive_numeric_order_gaps_and_repeat(self):
        self.fov(self.before, 'A-FOV')
        self.fov(self.before, 'A-FOV-13')
        for name in ['A-FOV', 'A-FOV-2', 'A-FOV-10', 'A-FOV-16', 'A-FOV-100']:
            self.fov(self.after, name, True)
        self.fov(self.after, 'B-FOV', True)
        self.fov(self.before, 'C-FOV-2')
        self.fov(self.after, 'C-FOV', True)
        before_snapshot = sorted(str(p) for p in self.root.rglob('*'))
        root, _, operations = self.plan()
        self.assertEqual(before_snapshot, sorted(str(p) for p in self.root.rglob('*')))
        oni.apply_plan(root, operations)
        for index, original in zip(range(14, 19), ['A-FOV', 'A-FOV-2', 'A-FOV-10', 'A-FOV-16', 'A-FOV-100']):
            name = f'A-FOV-{index}'
            nested = self.after / name / ('nested_' + name)
            self.assertEqual((nested / ('experiment_' + name + '_pos0.tif')).read_bytes(), b'raw\x00data')
            self.assertIn(original, (nested / 'metadata.json').read_text())
        self.assertTrue((self.after / 'B-FOV').exists())
        self.assertTrue((self.after / 'C-FOV-3').exists())
        self.assertTrue((self.before / 'A-FOV-13').exists())
        self.assertEqual(self.plan()[2], [])
        self.assertEqual(oni.main([str(self.before), str(self.after), '--apply']), 0)
        records = [json.loads(line) for line in (root / oni.JOURNAL).read_text().splitlines()]
        self.assertEqual(records[-1]['event'], 'complete')
        self.assertEqual(sum(r['event'] == 'renamed' for r in records), len(operations))

    def test_token_boundaries_and_symlinks(self):
        self.fov(self.before, 'A-FOV-2')
        bare = self.fov(self.after, 'A-FOV')
        numbered = self.fov(self.after, 'A-FOV-1')
        for path, names in [(bare, ['x_A-FOV.tif', 'x_A-FOV-1.tif', 'x_BA-FOV.tif']),
                            (numbered, ['x_A-FOV-1.tif', 'x_A-FOV-10.tif'])]:
            for name in names:
                (path / name).write_text('unchanged')
        (bare / 'link_A-FOV').symlink_to(self.before, target_is_directory=True)
        root, _, operations = self.plan()
        oni.apply_plan(root, operations)
        self.assertFalse((self.after / 'A-FOV-3/x_A-FOV.tif').exists())
        for name in ['x_A-FOV-3.tif', 'x_A-FOV-1.tif', 'x_BA-FOV.tif', 'link_A-FOV']:
            self.assertTrue((self.after / 'A-FOV-3' / name).exists())
        for name in ['x_A-FOV-4.tif', 'x_A-FOV-10.tif']:
            self.assertTrue((self.after / 'A-FOV-4' / name).exists())

    def test_collision_preflight(self):
        self.fov(self.before, 'A-FOV')
        path = self.fov(self.after, 'A-FOV')
        (path / 'A-FOV.tif').touch()
        (path / 'A-FOV-1.tif').touch()
        with self.assertRaisesRegex(ValueError, 'Conflicting destination'):
            self.plan()
        self.assertFalse((self.after / oni.JOURNAL).exists())
        self.assertTrue(path.exists())

    def test_invalid_roots_and_duplicates(self):
        for other in [self.before, self.root]:
            with self.assertRaises(ValueError):
                oni.build_plan(self.before, other)
        self.fov(self.after, 'A-FOV')
        self.fov(self.after, 'A-FOV-0')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.plan()

    def test_interruption_blocks_reapply(self):
        self.fov(self.before, 'A-FOV')
        self.fov(self.after, 'A-FOV', True)
        root, _, operations = self.plan()
        with patch.object(Path, 'rename', side_effect=OSError('simulated interruption')):
            with self.assertRaises(OSError):
                oni.apply_plan(root, operations)
        with self.assertRaisesRegex(ValueError, 'Incomplete journal'):
            oni.apply_plan(root, operations)
        self.assertTrue((self.after / 'A-FOV').exists())
        with self.assertRaisesRegex(ValueError, 'Incomplete journal'):
            oni.apply_plan(root, [])

    def test_spaces_in_condition_and_top_level_collision(self):
        self.fov(self.before, 'Condition A-FOV')
        self.fov(self.after, 'Condition A-FOV', True)
        (self.after / 'Condition A-FOV-1').write_text('occupied')
        with self.assertRaisesRegex(ValueError, 'Conflicting destination'):
            self.plan()
        (self.after / 'Condition A-FOV-1').unlink()
        root, _, operations = self.plan()
        oni.apply_plan(root, operations)
        self.assertTrue((self.after / 'Condition A-FOV-1').is_dir())


if __name__ == '__main__':
    unittest.main()
