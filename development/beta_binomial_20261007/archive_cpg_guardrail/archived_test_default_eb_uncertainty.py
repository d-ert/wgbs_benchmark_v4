"""Local CpG uncertainty must not inherit chromosome-wide residual df."""
import numpy as np
from epykit.dmc import _score_finalize


def _group(meth, depth=20):
    m = np.asarray(meth, dtype=float)
    n = np.full_like(m, depth)
    return n.sum(1), m.sum(1), (m*m/n).sum(1), np.full(len(m),m.shape[1],dtype=np.int32)


def _boundary_site(background_sites):
    background = np.tile([3,4,5,6,7], (background_sites,1))
    case = np.vstack([[2,3,4,5,6],background])
    control = np.vstack([[0,0,0,0,0],background])
    return _score_finalize(*_group(case), *_group(control), dispersion='eb')


def test_boundary_site_does_not_borrow_chromosome_residual_degrees_of_freedom():
    small = _boundary_site(100)
    large = _boundary_site(1000)
    # Adding identical background leaves the learned scale and prior weight
    # unchanged; it cannot turn a five-vs-five test into thousands of replicates.
    np.testing.assert_allclose(small[-1][0], large[-1][0])
    np.testing.assert_allclose(small[0][0], large[0][0])
    assert small[-1][0] < 100


def test_one_informative_group_contributes_only_its_residual_df():
    # Identical background gives the implementation's fixed prior weight=4.
    # Only the five case replicates contribute an estimable Pearson residual.
    result = _boundary_site(100)
    assert result[-1][0] == 4 + 4
    assert result[-1][1] == 8 + 4


def test_no_informative_group_uses_prior_uncertainty_only():
    background = np.tile([3,4,5,6,7], (100,1))
    case = np.vstack([[20]*5,background])
    control = np.vstack([[0]*5,background])
    result = _score_finalize(*_group(case), *_group(control), dispersion='eb')
    assert result[-1][0] == 4


def test_informative_boundary_group_retains_its_large_pearson_variance():
    background = np.tile([3,4,5,6,7], (100,1))
    case = np.vstack([[0,0,0,0,15],background])
    control = np.vstack([[0]*5,background])
    result = _score_finalize(*_group(case), *_group(control), dispersion='eb')
    # Q=(225/20 - 15**2/100)/(.15*.85), df=4, prior scale=1, weight=4.
    q = 9 / (.15*.85)
    expected_scale = (q + 4) / 8
    np.testing.assert_allclose(result[5][0], expected_scale)
    assert result[-1][0] == 8
    reverse = _score_finalize(*_group(control), *_group(case), dispersion='eb')
    for i in [0,5,6]:
        np.testing.assert_allclose(result[i], reverse[i])


def test_default_eb_maps_residual_scale_to_unequal_depth_contrast_variance():
    n = np.array([4,7,13,23,40], dtype=float)
    m1 = np.array([0,0,0,20,40], dtype=float)
    m2 = np.array([1,1,1,3,5], dtype=float)
    def summarize(m):
        return np.array([n.sum()]),np.array([m.sum()]),np.array([(m*m/n).sum()]),np.array([len(n)],dtype=np.int32)
    result = _score_finalize(*summarize(m1),*summarize(m2),dispersion='eb',
                             sn2_case=np.array([(n*n).sum()]),sn2_ctrl=np.array([(n*n).sum()]))
    # With one site, the unchanged prior fallback is scale=1, weight=4.
    def pearson(m):
        mu=m.sum()/n.sum()
        return np.sum((m-n*mu)**2/(n*mu*(1-mu)))
    residual_scale=(pearson(m1)+pearson(m2)+4)/(8+4)
    # Direct beta-binomial count variance for rho=1, normalized by the
    # binomial contrast variance; subtract one to obtain its rho coefficient.
    contrast_coeff=(2*np.sum(n*(n-1))/n.sum()**2)/(2/n.sum())
    residual_coeff=2*(n.sum()-len(n)-(np.sum(n*n)/n.sum()-1))/8
    expected=1+(residual_scale-1)*contrast_coeff/residual_coeff
    np.testing.assert_allclose(result[5][0],expected)
    reverse=_score_finalize(*summarize(m2),*summarize(m1),dispersion='eb',
                           sn2_case=np.array([(n*n).sum()]),sn2_ctrl=np.array([(n*n).sum()]))
    for i in [0,5,6]:np.testing.assert_allclose(result[i],reverse[i])


def test_equal_depth_count_geometry_preserves_the_eb_residual_scale():
    m1=np.array([[0,0,0,20,20]],dtype=float)
    m2=np.array([[1,1,1,3,5]],dtype=float)
    plain=_score_finalize(*_group(m1),*_group(m2),dispersion='eb')
    explicit=_score_finalize(*_group(m1),*_group(m2),dispersion='eb',
                            sn2_case=np.array([2000.]),sn2_ctrl=np.array([2000.]))
    for i in [0,5,6]:np.testing.assert_allclose(plain[i],explicit[i])
