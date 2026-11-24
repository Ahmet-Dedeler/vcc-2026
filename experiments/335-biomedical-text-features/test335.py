import numpy as np
import torch
from run335 import body,mean_pool


def test_identity():
    assert body({"A":"Gene Symbol A protein"},"A") == "protein"
    assert not body({"A":"Gene Symbol B protein"},"A")
    assert not body({"AA":"Gene Symbol AA protein"},"A")


def test_empty():
    assert not body({},"A")
    assert not body({"A":"Gene Symbol A"},"A")


def test_mean_pool_mask():
    hidden = torch.tensor([[[1.,3.],[3.,5.],[100.,100.]]])
    mask = torch.tensor([[1,1,0]])
    np.testing.assert_array_equal(mean_pool(hidden,mask).numpy(),[[2.,4.]])


def test_pool_padding_invariance():
    h = torch.tensor([[[1.,2.],[5.,6.]]])
    padded = torch.cat([h,torch.full((1,5,2),1000.)],dim=1)
    np.testing.assert_array_equal(mean_pool(h,torch.ones((1,2))).numpy(),
                                  mean_pool(padded,torch.tensor([[1,1,0,0,0,0,0]])).numpy())
