"""Regression tests: python -m unittest discover -s image_analysis/segmentation/tests."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    'sacd_nuclei', Path(__file__).resolve().parents[1] / 'cli_sacd_nuclei.py')
pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline)


def row(condition, mean, included=True):
    return dict(condition=condition, mean_spen_intensity=mean, included=included)


class ExpressionMatchingTests(unittest.TestCase):
    def test_arbitrary_conditions_and_inclusive_endpoints(self):
        rows = [row('SHA', 10), row('SHA', 20), row('FVP24h', 10),
                row('FVPrecover3h', 20), row('new-condition', 15),
                row('new-condition', 9), row('new-condition', 21),
                row('SHA', 1000, False), row('FVP24h', 15, False)]
        self.assertEqual(pipeline.match_expression(rows), [10, 20])
        self.assertEqual([r['expression_matched'] for r in rows],
                         [True, True, True, True, True, False, False, False, False])
        self.assertEqual(rows[5]['expression_exclusion_reason'], 'below_SHA_mean_spen_min')
        self.assertEqual(rows[6]['expression_exclusion_reason'], 'above_SHA_mean_spen_max')

    def test_configurable_reference_without_sha(self):
        rows = [row('control', 2), row('control', 4), row('treatment', 3), row('treatment', 5)]
        self.assertEqual(pipeline.match_expression(rows, 'control'), [2, 4])
        self.assertEqual(rows[-1]['expression_exclusion_reason'], 'above_control_mean_spen_max')
        self.assertTrue(rows[2]['expression_matched'])

    def test_missing_or_filtered_reference_reports_available_conditions(self):
        for rows in ([row('other', 5)], [row('SHA', 5, False), row('other', 6)]):
            with self.subTest(rows=rows):
                with self.assertRaisesRegex(ValueError, 'SHA.*Available conditions:.*other'):
                    pipeline.match_expression(rows)

    def test_preflight_missing_reference_before_output_creation(self):
        args = pipeline.parser().parse_args(['/does/not/exist', '--reference-condition', 'control'])
        with patch.object(pipeline, 'discover', return_value=(
                [{'condition': 'other', 'fov': 'other-FOV-1'}], ['DAPI', 'SPEN_JFX650'],
                {'unit': 'um', 'resolution': [1., 1.]})):
            with self.assertRaisesRegex(ValueError, "control.*Available conditions: other"):
                pipeline.preflight(args)

    def test_single_value_reference(self):
        rows = [row('SHA', 5), row('new', 5), row('new', 5.01)]
        self.assertEqual(pipeline.match_expression(rows), [5, 5])
        self.assertEqual([r['expression_matched'] for r in rows], [True, True, False])

    def test_panel_order_preserves_all_conditions_and_intensity_sort(self):
        names = ['dIDR', 'FVP24h', 'SHA', 'FVPrecover3h', 'dRRM']
        order = ['SHA', 'dRRM', 'dIDR']
        self.assertEqual(sorted(names, key=lambda c: pipeline.condition_order_key(c, order)),
                         ['SHA', 'dRRM', 'dIDR', 'FVP24h', 'FVPrecover3h'])
        rows = [dict(condition=c, mean_spen_intensity=v, expression_matched=True,
                     fov=f'{c}-FOV-1', mask_label=i)
                for c in names for i, v in enumerate([10, 20], 1)]
        result = pipeline.expression_sorted_rows(rows, order)
        self.assertEqual([r['condition'] for r in result[::2]],
                         ['SHA', 'dRRM', 'dIDR', 'FVP24h', 'FVPrecover3h'])
        self.assertEqual([r['mean_spen_intensity'] for r in result], [20, 10]*5)


if __name__ == '__main__':
    unittest.main()
