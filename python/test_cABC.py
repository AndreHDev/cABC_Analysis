import pytest
import numpy as np
from cABC_analysis import cABC_analysis

# ----------------- Tests for cABC_analysis -----------------

def test_identical_values():
    test1 = [1] * 13
    with pytest.warns(UserWarning, match="All data values are identical"):
        res = cABC_analysis(test1, plot_it=False)
    assert set(res['Aind']) == set(range(0, 13))
    assert res['Bind'] == []
    assert res['Cind'] == []

def test_character_input():
    test2 = ["10", 5, "x", 3, "0", 0]
    with pytest.warns(UserWarning, match="3 of 6 items are larger then 0"):
        res = cABC_analysis(test2, plot_it=False)
    assert set(res['Aind']) == {0}
    assert set(res['Bind']) == {1}
    assert set(res['Cind']) == {2, 3, 4, 5}

def test_a_lot_of_zeros():
    test3 = [0, 10, 0, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    with pytest.warns(UserWarning) as record:
        res = cABC_analysis(test3, plot_it=False)
    warning_messages = [str(w.message) for w in record]
    assert any("3 of 23 items are larger then 0" in msg for msg in warning_messages)
    assert any("duplicate value(s) spanning multiple classes" in msg for msg in warning_messages)

def test_AB_duplicate_values():
    test4 = [3, 1, 3, 1, 1, 3]
    with pytest.warns(UserWarning, match="duplicate value\\(s\\) spanning multiple classes"):
        res = cABC_analysis(test4, plot_it=False)
    assert set(res['Aind']) == {0, 2, 5}
    assert res['Bind'] == []
    assert set(res['Cind']) == {1, 3, 4}

# Diffrent from R version because of interpolation
def test_single_high_value_uniform_rest():
    test5 = [60] + [10] * 14
    with pytest.warns(UserWarning, match="duplicate value\\(s\\) spanning multiple classes"):
        res = cABC_analysis(test5, plot_it=False)
    assert res['Aind'] == [0]
    assert set(res['Bind']) == set(range(1, 15))
    assert res['Cind'] == []

def test_BC_duplicate_values():
    test6 = [10] * 3 + [5] * 3 + [1] * 3
    with pytest.warns(UserWarning, match="duplicate value\\(s\\) spanning multiple classes"):
        res = cABC_analysis(test6, plot_it=False)
    assert set(res['Aind']) == {0, 1, 2}
    assert set(res['Bind']) == {3, 4, 5}
    assert set(res['Cind']) == {6, 7, 8}

def test_linear_data():
    test7 = list(range(1, 11))
    # Expect no warnings
    res = cABC_analysis(test7, plot_it=False)
    # No specific assertions, just that it runs without warnings

def test_high_value_point():
    test8 = [100, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    res = cABC_analysis(test8, plot_it=False)
    assert res['Aind'] == [0]

def test_uniform_but_not_identical_data():
    np.random.seed(123)
    test9 = 100 + np.random.uniform(-0.0002, 0.0002, 100)
    # Expect no warnings
    res = cABC_analysis(test9, plot_it=False)

def test_uniform_but_not_identical_data_large_scale():
    np.random.seed(456)
    test10 = 100 + np.random.uniform(-0.0002, 0.0002, 10000)
    # Expect no warnings
    res = cABC_analysis(test10, plot_it=False)

def test_ABC_curve_under_uniform_at_beginning():
    test11 = [10, 10, 1, 1, 1]
    with pytest.warns(UserWarning, match="duplicate value\\(s\\) spanning multiple classes"):
        res = cABC_analysis(test11, plot_it=False)

def test_NULL_input():
    with pytest.raises(ValueError):
        cABC_analysis(None)

def test_single_number_input():
    with pytest.warns(UserWarning, match="Only one data point remains after filtering"):
        res = cABC_analysis(42, plot_it=False)
    assert res['Aind'] == [0]
    assert res['Bind'] == []
    assert res['Cind'] == []

def test_single_number_input_vector():
    with pytest.warns(UserWarning, match="Only one data point remains after filtering"):
        res = cABC_analysis([42], plot_it=False)
    assert res['Aind'] == [0]
    assert res['Bind'] == []
    assert res['Cind'] == []

def test_input_with_NA_values():
    test12 = [10, np.nan, 5, 7, np.nan, 3]
    with pytest.warns(UserWarning, match="4 of 6 items are larger then 0"):
        res = cABC_analysis(test12, plot_it=False)
    assert set(res['Aind']) == {0, 3}
    assert set(res['Bind']) == {2}
    assert set(res['Cind']) == {1, 4, 5}

def test_less_than_3_values():
    test13 = [10, 4]
    with pytest.warns(UserWarning, match="Extremly small dataset"):
        res = cABC_analysis(test13, plot_it=False)

def test_less_than_3_values_with_zeros():
    test14 = [0, 10, 0, 4, 0, 0]
    with pytest.warns(UserWarning, match="2 of 6 items are larger then 0"):
        res = cABC_analysis(test14, plot_it=False)

def test_duplicate_value_spanning_all_three_sets():
    test15 = [1] + [10] * 15
    with pytest.warns(UserWarning, match="duplicate value\\(s\\) spanning multiple classes"):
        res = cABC_analysis(test15, plot_it=False)