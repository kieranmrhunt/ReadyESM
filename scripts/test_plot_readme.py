"""Run with python -m unittest discover -s scripts -p test_plot_readme.py."""
import unittest
import numpy as np
from plot_readme import daily_mean, nearest


class PlotDataTests(unittest.TestCase):
    def test_daily_endpoints_and_missing_initial_sample(self):
        t = np.arange(0, 2.25, .25)
        v = np.array([np.nan, 1, 2, 3, 4, 5, 6, 7, 8])
        days, values = daily_mean(t, v, 2)
        np.testing.assert_array_equal(days, [1, 2])
        np.testing.assert_array_equal(values, [2.5, 6.5])

    def test_missing_sample_rejected(self):
        with self.assertRaises(ValueError):
            daily_mean(np.arange(0, 1.25, .25), np.array([0, 1, np.nan, 3, 4]), 1)

    def test_irregular_sample_rejected(self):
        with self.assertRaises(ValueError):
            daily_mean(np.array([0, .25, .5, 1]), np.arange(4), 1)

    def test_incomplete_day_rejected(self):
        with self.assertRaises(ValueError):
            daily_mean(np.arange(0, 1, .25), np.arange(4), 1)

    def test_periodic_spherical_neighbours(self):
        index = nearest(np.array([179, 0]), np.array([0, 0]),
                        np.array([[-179, 359]]), np.array([[0, 0]]))
        np.testing.assert_array_equal(index, [[0, 1]])

    def test_dry_cell_is_retained(self):
        index = nearest(np.array([0, 10]), np.array([0, 0]),
                        np.array([[1, 9]]), np.array([[0, 0]]))
        field = np.array([np.nan, 25])[index]
        self.assertTrue(np.isnan(field[0, 0]))
        self.assertEqual(field[0, 1], 25)


if __name__ == "__main__":
    unittest.main()
