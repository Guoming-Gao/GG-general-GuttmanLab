from types import SimpleNamespace

import numpy as np

from oni_spt.detection import _detail_value


def test_spotiflow_detail_value_accepts_scalar_and_single_channel_arrays():
    details = SimpleNamespace(
        prob=np.array([0.7, 0.8]),
        intens=np.array([[10.0], [12.0]]),
    )
    assert _detail_value(details, "prob", 1) == 0.8
    assert _detail_value(details, "intens", 1) == 12.0
    assert np.isnan(_detail_value(details, "missing", 0))
