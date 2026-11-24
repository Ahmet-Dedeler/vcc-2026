# C7 — saved-array observed-response readout

Validity: **VALID**. Four-cell Holm FWER 0.05.

C7a: observed median < 0 AND negative-tail excess > 0 AND Holm(max of two p).
C7b/c/d: positive observed rho in the declared direction AND Holm.

| Cell | Verdict | Raw p | Holm p | Estimable |
|---|---|---:|---:|---|
| C7a | PASS | 0.00049975 | 0.001999 | True |
| C7b | FAIL | 0.988006 | 1 | True |
| C7c | FAIL | 0.809595 | 1 | True |
| C7d | FAIL | 0.885057 | 1 | True |

95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI.

Undefined expression/rank cells remain in the family at p=1; no zero-signal imputation.
CPM premise failure voids every conclusion, regardless of nominal significance.

## Full readout, intervals, null summaries, provenance and limitations

```json
{
  "readout_version": "C7-fixed-axis-uniform-v2",
  "validity": {
    "cpm_gate_pass": true,
    "validity": "VALID",
    "shared_control_identical": true,
    "shared_control_interpretation": "retained",
    "cpm_rule": "all 510 source-axis CPM sums within 1 of 10000"
  },
  "cells": {
    "C7a": {
      "n": 100,
      "estimable": true,
      "reason": null,
      "observed": {
        "median": -0.026330770468529323,
        "tail_asym": 0.27999999999999997
      },
      "p_median_less": 0.0004997501249375312,
      "p_tail_greater": 0.0004997501249375312,
      "p_raw": 0.0004997501249375312,
      "direction_ok": true,
      "ci": {
        "median": {
          "lo": -0.05203105622228967,
          "hi": -0.009486476782116924,
          "n_finite": 20000,
          "n_draws": 20000,
          "scope": "95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI"
        },
        "tail_asym": {
          "lo": 0.14,
          "hi": 0.41000000000000003,
          "n_finite": 20000,
          "n_draws": 20000,
          "scope": "95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI"
        }
      },
      "null": {
        "median": {
          "n_total": 2000,
          "n_finite": 2000,
          "mean": 0.0060744184963344786,
          "sd": 0.006893038709532028,
          "min": -0.01749544196623171,
          "q025": -0.006975037407496198,
          "median": 0.006022575352479102,
          "q975": 0.01926114269316598,
          "max": 0.02975726413116027
        },
        "tail_asym": {
          "n_total": 2000,
          "n_finite": 2000,
          "mean": -0.06603,
          "sd": 0.05900372106909869,
          "min": -0.28,
          "q025": -0.18,
          "median": -0.060000000000000026,
          "q975": 0.05000000000000002,
          "max": 0.13999999999999999
        }
      },
      "p_holm": 0.001999000499750125,
      "reject": true,
      "verdict": "PASS",
      "statistical_verdict": "PASS"
    },
    "C7b": {
      "n": 97,
      "feature_missing": 3,
      "estimable": true,
      "reason": null,
      "observed": -0.1288265306122449,
      "p_raw": 0.9880059970014993,
      "direction_ok": false,
      "ci": {
        "lo": -0.31226149933625125,
        "hi": 0.060571387931268386,
        "n_finite": 20000,
        "n_draws": 20000,
        "scope": "95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI"
      },
      "null": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": 0.08024868504102674,
        "sd": 0.09657236659535348,
        "min": -0.2027009257311172,
        "q025": -0.10356781243425206,
        "median": 0.07965363980643804,
        "q975": 0.27094236534820104,
        "max": 0.3825347149168946
      },
      "p_holm": 1.0,
      "reject": false,
      "verdict": "FAIL",
      "statistical_verdict": "FAIL"
    },
    "C7c": {
      "n": 97,
      "feature_missing": 3,
      "estimable": true,
      "reason": null,
      "observed": 0.07605722701451714,
      "p_raw": 0.8095952023988006,
      "direction_ok": true,
      "ci": {
        "lo": -0.12188580905557224,
        "hi": 0.27024002766287974,
        "n_finite": 20000,
        "n_draws": 20000,
        "scope": "95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI"
      },
      "null": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": 0.1557192957079739,
        "sd": 0.086749171411323,
        "min": -0.17477119713864928,
        "q025": -0.0187812302756154,
        "median": 0.1603723963812329,
        "q975": 0.32173692930780545,
        "max": 0.4601830422890806
      },
      "p_holm": 1.0,
      "reject": false,
      "verdict": "FAIL",
      "statistical_verdict": "FAIL"
    },
    "C7d": {
      "n": 97,
      "feature_missing": 3,
      "estimable": true,
      "reason": null,
      "observed": 0.08069903219019567,
      "p_raw": 0.8850574712643678,
      "direction_ok": true,
      "ci": {
        "lo": -0.13476320519748272,
        "hi": 0.2874533271626702,
        "n_finite": 20000,
        "n_draws": 20000,
        "scope": "95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI"
      },
      "null": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": 0.18421035398695562,
        "sd": 0.08497401518766863,
        "min": -0.0975962549968441,
        "q025": 0.011566050389227856,
        "median": 0.18559330948874395,
        "q975": 0.3472089995792131,
        "max": 0.4572901325478645
      },
      "p_holm": 1.0,
      "reject": false,
      "verdict": "FAIL",
      "statistical_verdict": "FAIL"
    }
  },
  "design": {
    "n_tf": 100,
    "measured_axis": 15309,
    "scored_axis": 15211,
    "dropped_own_columns": 98,
    "axis_sha256": "4f07baf5d234c31c059a8471577c0ba5a85640586dfd6f392cfc75ee4f98b428",
    "cohort_sha256": "fca12254412534255a618dc07aa9f2da3fbdcef05b4c00de6fe900efed0a07ff",
    "seed": 20260907,
    "n_null": 2000,
    "n_boot": 20000,
    "fixed_points": 2031,
    "permutation_scope": "uniform whole-cohort OE matching; fixed points allowed; KD/features fixed",
    "ci_scope": "95% percentile TF bootstrap, conditional on frozen centred rows and each cell's fixed feature-complete cohort; not between-study CI"
  },
  "observed": {
    "n": 100,
    "n_total": 100,
    "n_undefined": 0,
    "median": -0.026330770468529323,
    "mean": -0.05079188848871492,
    "share_neg": 0.64,
    "share_below": 0.41,
    "share_above": 0.13,
    "n_below": 41,
    "n_above": 13,
    "tail_asym": 0.27999999999999997
  },
  "secondaries": {
    "uncentered_fixed_axis": {
      "observed": {
        "n": 100,
        "n_total": 100,
        "n_undefined": 0,
        "median": -0.020291057855196717,
        "mean": -0.04454366410183911,
        "share_neg": 0.68,
        "share_below": 0.39,
        "share_above": 0.12,
        "n_below": 39,
        "n_above": 12,
        "tail_asym": 0.27
      },
      "null_median": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": 0.0007196711813570108,
        "sd": 0.004072072818838295,
        "min": -0.01232790077644599,
        "q025": -0.007127197631327703,
        "median": 0.000605914094470096,
        "q975": 0.008644120515562162,
        "max": 0.014100152971818802
      },
      "null_tail": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": -0.03882000000000001,
        "sd": 0.0438315822210424,
        "min": -0.18,
        "q025": -0.12000000000000001,
        "median": -0.039999999999999994,
        "q975": 0.04999999999999999,
        "max": 0.11000000000000001
      },
      "verdict": "DESCRIPTIVE ONLY"
    },
    "bankwide_centered_fixed_axis": {
      "observed": {
        "n": 100,
        "n_total": 100,
        "n_undefined": 0,
        "median": -0.03705186954327948,
        "mean": -0.04543062867277966,
        "share_neg": 0.6,
        "share_below": 0.44,
        "share_above": 0.15,
        "n_below": 44,
        "n_above": 15,
        "tail_asym": 0.29000000000000004
      },
      "null_median": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": 0.011236504538015153,
        "sd": 0.006588202624648157,
        "min": -0.010772580692058489,
        "q025": -0.0014183963509934853,
        "median": 0.011191363927064081,
        "q975": 0.024270656934737924,
        "max": 0.034236885987083396
      },
      "null_tail": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": -0.102625,
        "sd": 0.05364568365674912,
        "min": -0.27,
        "q025": -0.21,
        "median": -0.1,
        "q975": 0.0,
        "max": 0.07000000000000003
      },
      "verdict": "DESCRIPTIVE ONLY"
    },
    "original_pairwise_own_exclusion": {
      "observed": {
        "n": 100,
        "n_total": 100,
        "n_undefined": 0,
        "median": -0.025539753574653698,
        "mean": -0.05098884866370968,
        "share_neg": 0.65,
        "share_below": 0.41,
        "share_above": 0.13,
        "n_below": 41,
        "n_above": 13,
        "tail_asym": 0.27999999999999997
      },
      "null_median": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": 0.006096898074232967,
        "sd": 0.006910235631699952,
        "min": -0.016986477894123955,
        "q025": -0.006827345761179096,
        "median": 0.006062682489903944,
        "q975": 0.01962672239805858,
        "max": 0.02962073114210387
      },
      "null_tail": {
        "n_total": 2000,
        "n_finite": 2000,
        "mean": -0.06698,
        "sd": 0.05932267357427513,
        "min": -0.28,
        "q025": -0.18,
        "median": -0.06999999999999998,
        "q975": 0.05000000000000002,
        "max": 0.13999999999999999
      },
      "verdict": "DESCRIPTIVE ONLY"
    },
    "own_transcript": {
      "n": 98,
      "oe_median": 2.25782310962677,
      "kd_median": -0.6734746396541595,
      "oe_share_up": 0.8571428571428571,
      "kd_share_down": 1.0,
      "S_on_own_up": {
        "n": 84,
        "n_total": 84,
        "n_undefined": 0,
        "median": -0.023069282366390498,
        "mean": -0.04395008647175332,
        "share_neg": 0.6190476190476191,
        "share_below": 0.40476190476190477,
        "share_above": 0.14285714285714285,
        "n_below": 34,
        "n_above": 12,
        "tail_asym": 0.2619047619047619
      },
      "db_joung_cohort_overlap_n": 1,
      "S_on_db_joung_cohort_overlap": {
        "n": 1,
        "n_total": 1,
        "n_undefined": 0,
        "median": -0.04511697841159331,
        "mean": -0.04511697841159331,
        "share_neg": 1.0,
        "share_below": 0.0,
        "share_above": 0.0,
        "n_below": 0,
        "n_above": 0,
        "tail_asym": 0.0
      }
    },
    "302": {
      "ALL58": {
        "n": 58,
        "rho": 0.4261281491279338
      },
      "primary45": {
        "n": 45,
        "rho": 0.5023715415019763
      },
      "scope": "saved SYM302 S_g versus C7 S; descriptive; C6 processing VOID; no reliability plug-in"
    },
    "baselines": {
      "C_excl": {
        "minusS": {
          "n": 97,
          "rho": -0.1288265306122449
        },
        "absS": {
          "n": 97,
          "rho": -0.12071323374710709
        }
      },
      "C_held": {
        "minusS": {
          "n": 97,
          "rho": -0.11266568483063329
        },
        "absS": {
          "n": 97,
          "rho": -0.09470334525562803
        }
      },
      "gen_pc20": {
        "minusS": {
          "n": 97,
          "rho": 0.07605722701451714
        },
        "absS": {
          "n": 97,
          "rho": 0.08069903219019567
        }
      },
      "breadth": {
        "minusS": {
          "n": 97,
          "rho": 0.017645227011345754
        },
        "absS": {
          "n": 97,
          "rho": -0.1459633728250964
        }
      },
      "mag290": {
        "minusS": {
          "n": 97,
          "rho": 0.009914790658531453
        },
        "absS": {
          "n": 97,
          "rho": -0.15801861981906165
        }
      },
      "r_bar": {
        "minusS": {
          "n": 97,
          "rho": -0.14726225541763097
        },
        "absS": {
          "n": 97,
          "rho": -0.1171496949295182
        }
      },
      "KD_magnitude": {
        "minusS": {
          "n": 100,
          "rho": -0.011329132913291328
        },
        "absS": {
          "n": 100,
          "rho": -0.1635763576357636
        }
      },
      "OE_magnitude": {
        "minusS": {
          "n": 100,
          "rho": 0.3742214221422142
        },
        "absS": {
          "n": 100,
          "rho": 0.2593939393939394
        }
      }
    },
    "partial_rho_given_magnitudes": {
      "C7b": -0.07797152708624994,
      "C7c": 0.1380321984239146,
      "C7d": 0.16952487254005932
    }
  },
  "no_signal_rows": {
    "OE": [],
    "KD": []
  },
  "feature_provenance": {
    "excluded_source": "gse295214_h1",
    "n_pairs": 338,
    "pairs": [
      [
        "altair_a549",
        "dld1_lowmoi"
      ],
      [
        "altair_a549",
        "frangieh_melanoma"
      ],
      [
        "altair_a549",
        "gasperini_k562"
      ],
      [
        "altair_a549",
        "h1_2025"
      ],
      [
        "altair_a549",
        "igvf_cardiomyocyte"
      ],
      [
        "altair_a549",
        "igvf_telohaec"
      ],
      [
        "altair_a549",
        "jiang_BXPC3"
      ],
      [
        "altair_a549",
        "jiang_HAP1"
      ],
      [
        "altair_a549",
        "jiang_HT29"
      ],
      [
        "altair_a549",
        "jiang_K562"
      ],
      [
        "altair_a549",
        "jiang_MCF7"
      ],
      [
        "altair_a549",
        "k562_essential"
      ],
      [
        "altair_a549",
        "k562_gwps"
      ],
      [
        "altair_a549",
        "nadig_hepg2"
      ],
      [
        "altair_a549",
        "nadig_jurkat"
      ],
      [
        "altair_a549",
        "orion_hct116"
      ],
      [
        "altair_a549",
        "orion_hek293t"
      ],
      [
        "altair_a549",
        "perturbmatch"
      ],
      [
        "altair_a549",
        "rpe1"
      ],
      [
        "altair_a549",
        "sunshine_calu3"
      ],
      [
        "altair_a549",
        "tian2019_ipsc"
      ],
      [
        "altair_a549",
        "tian_ipsc_neuron"
      ],
      [
        "altair_a549",
        "wessels_thp1"
      ],
      [
        "altair_a549",
        "xucao_hek293"
      ],
      [
        "dld1_lowmoi",
        "frangieh_melanoma"
      ],
      [
        "dld1_lowmoi",
        "gasperini_k562"
      ],
      [
        "dld1_lowmoi",
        "h1_2025"
      ],
      [
        "dld1_lowmoi",
        "igvf_cardiomyocyte"
      ],
      [
        "dld1_lowmoi",
        "igvf_telohaec"
      ],
      [
        "dld1_lowmoi",
        "jiang_A549"
      ],
      [
        "dld1_lowmoi",
        "jiang_BXPC3"
      ],
      [
        "dld1_lowmoi",
        "jiang_HAP1"
      ],
      [
        "dld1_lowmoi",
        "jiang_HT29"
      ],
      [
        "dld1_lowmoi",
        "jiang_K562"
      ],
      [
        "dld1_lowmoi",
        "jiang_MCF7"
      ],
      [
        "dld1_lowmoi",
        "k562_essential"
      ],
      [
        "dld1_lowmoi",
        "k562_gwps"
      ],
      [
        "dld1_lowmoi",
        "myllia_a549"
      ],
      [
        "dld1_lowmoi",
        "nadig_hepg2"
      ],
      [
        "dld1_lowmoi",
        "nadig_jurkat"
      ],
      [
        "dld1_lowmoi",
        "orion_hct116"
      ],
      [
        "dld1_lowmoi",
        "orion_hek293t"
      ],
      [
        "dld1_lowmoi",
        "perturbmatch"
      ],
      [
        "dld1_lowmoi",
        "rpe1"
      ],
      [
        "dld1_lowmoi",
        "sunshine_calu3"
      ],
      [
        "dld1_lowmoi",
        "tian2019_ipsc"
      ],
      [
        "dld1_lowmoi",
        "tian_ipsc_neuron"
      ],
      [
        "dld1_lowmoi",
        "wessels_thp1"
      ],
      [
        "dld1_lowmoi",
        "xucao_hek293"
      ],
      [
        "frangieh_melanoma",
        "gasperini_k562"
      ],
      [
        "frangieh_melanoma",
        "h1_2025"
      ],
      [
        "frangieh_melanoma",
        "igvf_cardiomyocyte"
      ],
      [
        "frangieh_melanoma",
        "igvf_telohaec"
      ],
      [
        "frangieh_melanoma",
        "jiang_A549"
      ],
      [
        "frangieh_melanoma",
        "jiang_BXPC3"
      ],
      [
        "frangieh_melanoma",
        "jiang_HAP1"
      ],
      [
        "frangieh_melanoma",
        "jiang_HT29"
      ],
      [
        "frangieh_melanoma",
        "jiang_K562"
      ],
      [
        "frangieh_melanoma",
        "jiang_MCF7"
      ],
      [
        "frangieh_melanoma",
        "k562_essential"
      ],
      [
        "frangieh_melanoma",
        "k562_gwps"
      ],
      [
        "frangieh_melanoma",
        "myllia_a549"
      ],
      [
        "frangieh_melanoma",
        "nadig_hepg2"
      ],
      [
        "frangieh_melanoma",
        "nadig_jurkat"
      ],
      [
        "frangieh_melanoma",
        "orion_hct116"
      ],
      [
        "frangieh_melanoma",
        "orion_hek293t"
      ],
      [
        "frangieh_melanoma",
        "perturbmatch"
      ],
      [
        "frangieh_melanoma",
        "rpe1"
      ],
      [
        "frangieh_melanoma",
        "sunshine_calu3"
      ],
      [
        "frangieh_melanoma",
        "tian2019_ipsc"
      ],
      [
        "frangieh_melanoma",
        "tian_ipsc_neuron"
      ],
      [
        "frangieh_melanoma",
        "wessels_thp1"
      ],
      [
        "frangieh_melanoma",
        "xucao_hek293"
      ],
      [
        "gasperini_k562",
        "h1_2025"
      ],
      [
        "gasperini_k562",
        "igvf_cardiomyocyte"
      ],
      [
        "gasperini_k562",
        "igvf_telohaec"
      ],
      [
        "gasperini_k562",
        "jiang_A549"
      ],
      [
        "gasperini_k562",
        "jiang_BXPC3"
      ],
      [
        "gasperini_k562",
        "jiang_HAP1"
      ],
      [
        "gasperini_k562",
        "jiang_HT29"
      ],
      [
        "gasperini_k562",
        "jiang_MCF7"
      ],
      [
        "gasperini_k562",
        "myllia_a549"
      ],
      [
        "gasperini_k562",
        "nadig_hepg2"
      ],
      [
        "gasperini_k562",
        "nadig_jurkat"
      ],
      [
        "gasperini_k562",
        "orion_hct116"
      ],
      [
        "gasperini_k562",
        "orion_hek293t"
      ],
      [
        "gasperini_k562",
        "rpe1"
      ],
      [
        "gasperini_k562",
        "sunshine_calu3"
      ],
      [
        "gasperini_k562",
        "tian2019_ipsc"
      ],
      [
        "gasperini_k562",
        "tian_ipsc_neuron"
      ],
      [
        "gasperini_k562",
        "wessels_thp1"
      ],
      [
        "gasperini_k562",
        "xucao_hek293"
      ],
      [
        "h1_2025",
        "igvf_cardiomyocyte"
      ],
      [
        "h1_2025",
        "igvf_telohaec"
      ],
      [
        "h1_2025",
        "jiang_A549"
      ],
      [
        "h1_2025",
        "jiang_BXPC3"
      ],
      [
        "h1_2025",
        "jiang_HAP1"
      ],
      [
        "h1_2025",
        "jiang_HT29"
      ],
      [
        "h1_2025",
        "jiang_K562"
      ],
      [
        "h1_2025",
        "jiang_MCF7"
      ],
      [
        "h1_2025",
        "k562_essential"
      ],
      [
        "h1_2025",
        "k562_gwps"
      ],
      [
        "h1_2025",
        "myllia_a549"
      ],
      [
        "h1_2025",
        "nadig_hepg2"
      ],
      [
        "h1_2025",
        "nadig_jurkat"
      ],
      [
        "h1_2025",
        "orion_hct116"
      ],
      [
        "h1_2025",
        "orion_hek293t"
      ],
      [
        "h1_2025",
        "perturbmatch"
      ],
      [
        "h1_2025",
        "rpe1"
      ],
      [
        "h1_2025",
        "sunshine_calu3"
      ],
      [
        "h1_2025",
        "tian2019_ipsc"
      ],
      [
        "h1_2025",
        "tian_ipsc_neuron"
      ],
      [
        "h1_2025",
        "wessels_thp1"
      ],
      [
        "h1_2025",
        "xucao_hek293"
      ],
      [
        "igvf_cardiomyocyte",
        "igvf_telohaec"
      ],
      [
        "igvf_cardiomyocyte",
        "jiang_A549"
      ],
      [
        "igvf_cardiomyocyte",
        "jiang_BXPC3"
      ],
      [
        "igvf_cardiomyocyte",
        "jiang_HAP1"
      ],
      [
        "igvf_cardiomyocyte",
        "jiang_HT29"
      ],
      [
        "igvf_cardiomyocyte",
        "jiang_K562"
      ],
      [
        "igvf_cardiomyocyte",
        "jiang_MCF7"
      ],
      [
        "igvf_cardiomyocyte",
        "k562_essential"
      ],
      [
        "igvf_cardiomyocyte",
        "k562_gwps"
      ],
      [
        "igvf_cardiomyocyte",
        "myllia_a549"
      ],
      [
        "igvf_cardiomyocyte",
        "nadig_hepg2"
      ],
      [
        "igvf_cardiomyocyte",
        "nadig_jurkat"
      ],
      [
        "igvf_cardiomyocyte",
        "orion_hct116"
      ],
      [
        "igvf_cardiomyocyte",
        "orion_hek293t"
      ],
      [
        "igvf_cardiomyocyte",
        "perturbmatch"
      ],
      [
        "igvf_cardiomyocyte",
        "rpe1"
      ],
      [
        "igvf_cardiomyocyte",
        "sunshine_calu3"
      ],
      [
        "igvf_cardiomyocyte",
        "tian2019_ipsc"
      ],
      [
        "igvf_cardiomyocyte",
        "tian_ipsc_neuron"
      ],
      [
        "igvf_cardiomyocyte",
        "wessels_thp1"
      ],
      [
        "igvf_cardiomyocyte",
        "xucao_hek293"
      ],
      [
        "igvf_telohaec",
        "jiang_A549"
      ],
      [
        "igvf_telohaec",
        "jiang_BXPC3"
      ],
      [
        "igvf_telohaec",
        "jiang_HAP1"
      ],
      [
        "igvf_telohaec",
        "jiang_HT29"
      ],
      [
        "igvf_telohaec",
        "jiang_K562"
      ],
      [
        "igvf_telohaec",
        "jiang_MCF7"
      ],
      [
        "igvf_telohaec",
        "k562_essential"
      ],
      [
        "igvf_telohaec",
        "k562_gwps"
      ],
      [
        "igvf_telohaec",
        "myllia_a549"
      ],
      [
        "igvf_telohaec",
        "nadig_hepg2"
      ],
      [
        "igvf_telohaec",
        "nadig_jurkat"
      ],
      [
        "igvf_telohaec",
        "orion_hct116"
      ],
      [
        "igvf_telohaec",
        "orion_hek293t"
      ],
      [
        "igvf_telohaec",
        "perturbmatch"
      ],
      [
        "igvf_telohaec",
        "rpe1"
      ],
      [
        "igvf_telohaec",
        "sunshine_calu3"
      ],
      [
        "igvf_telohaec",
        "tian2019_ipsc"
      ],
      [
        "igvf_telohaec",
        "tian_ipsc_neuron"
      ],
      [
        "igvf_telohaec",
        "wessels_thp1"
      ],
      [
        "igvf_telohaec",
        "xucao_hek293"
      ],
      [
        "jiang_A549",
        "jiang_BXPC3"
      ],
      [
        "jiang_A549",
        "jiang_HAP1"
      ],
      [
        "jiang_A549",
        "jiang_HT29"
      ],
      [
        "jiang_A549",
        "jiang_K562"
      ],
      [
        "jiang_A549",
        "jiang_MCF7"
      ],
      [
        "jiang_A549",
        "k562_essential"
      ],
      [
        "jiang_A549",
        "k562_gwps"
      ],
      [
        "jiang_A549",
        "nadig_hepg2"
      ],
      [
        "jiang_A549",
        "nadig_jurkat"
      ],
      [
        "jiang_A549",
        "orion_hct116"
      ],
      [
        "jiang_A549",
        "orion_hek293t"
      ],
      [
        "jiang_A549",
        "perturbmatch"
      ],
      [
        "jiang_A549",
        "rpe1"
      ],
      [
        "jiang_A549",
        "sunshine_calu3"
      ],
      [
        "jiang_A549",
        "tian2019_ipsc"
      ],
      [
        "jiang_A549",
        "tian_ipsc_neuron"
      ],
      [
        "jiang_A549",
        "wessels_thp1"
      ],
      [
        "jiang_A549",
        "xucao_hek293"
      ],
      [
        "jiang_BXPC3",
        "jiang_HAP1"
      ],
      [
        "jiang_BXPC3",
        "jiang_HT29"
      ],
      [
        "jiang_BXPC3",
        "jiang_K562"
      ],
      [
        "jiang_BXPC3",
        "jiang_MCF7"
      ],
      [
        "jiang_BXPC3",
        "k562_essential"
      ],
      [
        "jiang_BXPC3",
        "k562_gwps"
      ],
      [
        "jiang_BXPC3",
        "myllia_a549"
      ],
      [
        "jiang_BXPC3",
        "nadig_hepg2"
      ],
      [
        "jiang_BXPC3",
        "nadig_jurkat"
      ],
      [
        "jiang_BXPC3",
        "orion_hct116"
      ],
      [
        "jiang_BXPC3",
        "orion_hek293t"
      ],
      [
        "jiang_BXPC3",
        "perturbmatch"
      ],
      [
        "jiang_BXPC3",
        "rpe1"
      ],
      [
        "jiang_BXPC3",
        "sunshine_calu3"
      ],
      [
        "jiang_BXPC3",
        "tian2019_ipsc"
      ],
      [
        "jiang_BXPC3",
        "tian_ipsc_neuron"
      ],
      [
        "jiang_BXPC3",
        "wessels_thp1"
      ],
      [
        "jiang_BXPC3",
        "xucao_hek293"
      ],
      [
        "jiang_HAP1",
        "jiang_HT29"
      ],
      [
        "jiang_HAP1",
        "jiang_K562"
      ],
      [
        "jiang_HAP1",
        "jiang_MCF7"
      ],
      [
        "jiang_HAP1",
        "k562_essential"
      ],
      [
        "jiang_HAP1",
        "k562_gwps"
      ],
      [
        "jiang_HAP1",
        "myllia_a549"
      ],
      [
        "jiang_HAP1",
        "nadig_hepg2"
      ],
      [
        "jiang_HAP1",
        "nadig_jurkat"
      ],
      [
        "jiang_HAP1",
        "orion_hct116"
      ],
      [
        "jiang_HAP1",
        "orion_hek293t"
      ],
      [
        "jiang_HAP1",
        "perturbmatch"
      ],
      [
        "jiang_HAP1",
        "rpe1"
      ],
      [
        "jiang_HAP1",
        "sunshine_calu3"
      ],
      [
        "jiang_HAP1",
        "tian2019_ipsc"
      ],
      [
        "jiang_HAP1",
        "tian_ipsc_neuron"
      ],
      [
        "jiang_HAP1",
        "wessels_thp1"
      ],
      [
        "jiang_HAP1",
        "xucao_hek293"
      ],
      [
        "jiang_HT29",
        "jiang_K562"
      ],
      [
        "jiang_HT29",
        "jiang_MCF7"
      ],
      [
        "jiang_HT29",
        "k562_essential"
      ],
      [
        "jiang_HT29",
        "k562_gwps"
      ],
      [
        "jiang_HT29",
        "myllia_a549"
      ],
      [
        "jiang_HT29",
        "nadig_hepg2"
      ],
      [
        "jiang_HT29",
        "nadig_jurkat"
      ],
      [
        "jiang_HT29",
        "orion_hct116"
      ],
      [
        "jiang_HT29",
        "orion_hek293t"
      ],
      [
        "jiang_HT29",
        "perturbmatch"
      ],
      [
        "jiang_HT29",
        "rpe1"
      ],
      [
        "jiang_HT29",
        "sunshine_calu3"
      ],
      [
        "jiang_HT29",
        "tian2019_ipsc"
      ],
      [
        "jiang_HT29",
        "tian_ipsc_neuron"
      ],
      [
        "jiang_HT29",
        "wessels_thp1"
      ],
      [
        "jiang_HT29",
        "xucao_hek293"
      ],
      [
        "jiang_K562",
        "jiang_MCF7"
      ],
      [
        "jiang_K562",
        "myllia_a549"
      ],
      [
        "jiang_K562",
        "nadig_hepg2"
      ],
      [
        "jiang_K562",
        "nadig_jurkat"
      ],
      [
        "jiang_K562",
        "orion_hct116"
      ],
      [
        "jiang_K562",
        "orion_hek293t"
      ],
      [
        "jiang_K562",
        "rpe1"
      ],
      [
        "jiang_K562",
        "sunshine_calu3"
      ],
      [
        "jiang_K562",
        "tian2019_ipsc"
      ],
      [
        "jiang_K562",
        "tian_ipsc_neuron"
      ],
      [
        "jiang_K562",
        "wessels_thp1"
      ],
      [
        "jiang_K562",
        "xucao_hek293"
      ],
      [
        "jiang_MCF7",
        "k562_essential"
      ],
      [
        "jiang_MCF7",
        "k562_gwps"
      ],
      [
        "jiang_MCF7",
        "myllia_a549"
      ],
      [
        "jiang_MCF7",
        "nadig_hepg2"
      ],
      [
        "jiang_MCF7",
        "nadig_jurkat"
      ],
      [
        "jiang_MCF7",
        "orion_hct116"
      ],
      [
        "jiang_MCF7",
        "orion_hek293t"
      ],
      [
        "jiang_MCF7",
        "perturbmatch"
      ],
      [
        "jiang_MCF7",
        "rpe1"
      ],
      [
        "jiang_MCF7",
        "sunshine_calu3"
      ],
      [
        "jiang_MCF7",
        "tian2019_ipsc"
      ],
      [
        "jiang_MCF7",
        "tian_ipsc_neuron"
      ],
      [
        "jiang_MCF7",
        "wessels_thp1"
      ],
      [
        "jiang_MCF7",
        "xucao_hek293"
      ],
      [
        "k562_essential",
        "myllia_a549"
      ],
      [
        "k562_essential",
        "nadig_hepg2"
      ],
      [
        "k562_essential",
        "nadig_jurkat"
      ],
      [
        "k562_essential",
        "orion_hct116"
      ],
      [
        "k562_essential",
        "orion_hek293t"
      ],
      [
        "k562_essential",
        "rpe1"
      ],
      [
        "k562_essential",
        "sunshine_calu3"
      ],
      [
        "k562_essential",
        "tian2019_ipsc"
      ],
      [
        "k562_essential",
        "tian_ipsc_neuron"
      ],
      [
        "k562_essential",
        "wessels_thp1"
      ],
      [
        "k562_essential",
        "xucao_hek293"
      ],
      [
        "k562_gwps",
        "myllia_a549"
      ],
      [
        "k562_gwps",
        "nadig_hepg2"
      ],
      [
        "k562_gwps",
        "nadig_jurkat"
      ],
      [
        "k562_gwps",
        "orion_hct116"
      ],
      [
        "k562_gwps",
        "orion_hek293t"
      ],
      [
        "k562_gwps",
        "rpe1"
      ],
      [
        "k562_gwps",
        "sunshine_calu3"
      ],
      [
        "k562_gwps",
        "tian2019_ipsc"
      ],
      [
        "k562_gwps",
        "tian_ipsc_neuron"
      ],
      [
        "k562_gwps",
        "wessels_thp1"
      ],
      [
        "k562_gwps",
        "xucao_hek293"
      ],
      [
        "myllia_a549",
        "nadig_hepg2"
      ],
      [
        "myllia_a549",
        "nadig_jurkat"
      ],
      [
        "myllia_a549",
        "orion_hct116"
      ],
      [
        "myllia_a549",
        "orion_hek293t"
      ],
      [
        "myllia_a549",
        "perturbmatch"
      ],
      [
        "myllia_a549",
        "rpe1"
      ],
      [
        "myllia_a549",
        "sunshine_calu3"
      ],
      [
        "myllia_a549",
        "tian2019_ipsc"
      ],
      [
        "myllia_a549",
        "tian_ipsc_neuron"
      ],
      [
        "myllia_a549",
        "wessels_thp1"
      ],
      [
        "myllia_a549",
        "xucao_hek293"
      ],
      [
        "nadig_hepg2",
        "nadig_jurkat"
      ],
      [
        "nadig_hepg2",
        "orion_hct116"
      ],
      [
        "nadig_hepg2",
        "orion_hek293t"
      ],
      [
        "nadig_hepg2",
        "perturbmatch"
      ],
      [
        "nadig_hepg2",
        "rpe1"
      ],
      [
        "nadig_hepg2",
        "sunshine_calu3"
      ],
      [
        "nadig_hepg2",
        "tian2019_ipsc"
      ],
      [
        "nadig_hepg2",
        "tian_ipsc_neuron"
      ],
      [
        "nadig_hepg2",
        "wessels_thp1"
      ],
      [
        "nadig_hepg2",
        "xucao_hek293"
      ],
      [
        "nadig_jurkat",
        "orion_hct116"
      ],
      [
        "nadig_jurkat",
        "orion_hek293t"
      ],
      [
        "nadig_jurkat",
        "perturbmatch"
      ],
      [
        "nadig_jurkat",
        "rpe1"
      ],
      [
        "nadig_jurkat",
        "sunshine_calu3"
      ],
      [
        "nadig_jurkat",
        "tian2019_ipsc"
      ],
      [
        "nadig_jurkat",
        "tian_ipsc_neuron"
      ],
      [
        "nadig_jurkat",
        "wessels_thp1"
      ],
      [
        "nadig_jurkat",
        "xucao_hek293"
      ],
      [
        "orion_hct116",
        "orion_hek293t"
      ],
      [
        "orion_hct116",
        "perturbmatch"
      ],
      [
        "orion_hct116",
        "rpe1"
      ],
      [
        "orion_hct116",
        "sunshine_calu3"
      ],
      [
        "orion_hct116",
        "tian2019_ipsc"
      ],
      [
        "orion_hct116",
        "tian_ipsc_neuron"
      ],
      [
        "orion_hct116",
        "wessels_thp1"
      ],
      [
        "orion_hct116",
        "xucao_hek293"
      ],
      [
        "orion_hek293t",
        "perturbmatch"
      ],
      [
        "orion_hek293t",
        "rpe1"
      ],
      [
        "orion_hek293t",
        "sunshine_calu3"
      ],
      [
        "orion_hek293t",
        "tian2019_ipsc"
      ],
      [
        "orion_hek293t",
        "tian_ipsc_neuron"
      ],
      [
        "orion_hek293t",
        "wessels_thp1"
      ],
      [
        "orion_hek293t",
        "xucao_hek293"
      ],
      [
        "perturbmatch",
        "rpe1"
      ],
      [
        "perturbmatch",
        "sunshine_calu3"
      ],
      [
        "perturbmatch",
        "tian2019_ipsc"
      ],
      [
        "perturbmatch",
        "tian_ipsc_neuron"
      ],
      [
        "perturbmatch",
        "wessels_thp1"
      ],
      [
        "perturbmatch",
        "xucao_hek293"
      ],
      [
        "rpe1",
        "sunshine_calu3"
      ],
      [
        "rpe1",
        "tian2019_ipsc"
      ],
      [
        "rpe1",
        "tian_ipsc_neuron"
      ],
      [
        "rpe1",
        "wessels_thp1"
      ],
      [
        "rpe1",
        "xucao_hek293"
      ],
      [
        "sunshine_calu3",
        "tian2019_ipsc"
      ],
      [
        "sunshine_calu3",
        "tian_ipsc_neuron"
      ],
      [
        "sunshine_calu3",
        "wessels_thp1"
      ],
      [
        "sunshine_calu3",
        "xucao_hek293"
      ],
      [
        "tian2019_ipsc",
        "tian_ipsc_neuron"
      ],
      [
        "tian2019_ipsc",
        "wessels_thp1"
      ],
      [
        "tian2019_ipsc",
        "xucao_hek293"
      ],
      [
        "tian_ipsc_neuron",
        "wessels_thp1"
      ],
      [
        "tian_ipsc_neuron",
        "xucao_hek293"
      ],
      [
        "wessels_thp1",
        "xucao_hek293"
      ]
    ],
    "n_pairs_per_tf": [
      33,
      9,
      63,
      25,
      3,
      33,
      33,
      42,
      63,
      42,
      27,
      9,
      9,
      27,
      10,
      9,
      20,
      63,
      52,
      20,
      33,
      12,
      33,
      33,
      87,
      52,
      9,
      113,
      63,
      33,
      33,
      44,
      9,
      14,
      9,
      75,
      9,
      5,
      42,
      20,
      20,
      25,
      126,
      14,
      33,
      33,
      25,
      33,
      14,
      0,
      20,
      10,
      9,
      42,
      52,
      25,
      25,
      35,
      14,
      21,
      5,
      9,
      63,
      25,
      87,
      35,
      63,
      51,
      25,
      0,
      65,
      42,
      146,
      25,
      88,
      9,
      20,
      20,
      99,
      14,
      14,
      20,
      9,
      0,
      6,
      36,
      63,
      9,
      42,
      42,
      18,
      20,
      25,
      9,
      33,
      52,
      20,
      25,
      9,
      21
    ]
  },
  "limitations": [
    "Cross-study/cross-line observed alignment, conditional on target-label exchangeability; not a causal null.",
    "No OE replicate structure: no disattenuation, latent-equivalence bound or reliability plug-in.",
    "HGNC family-block bootstrap NOT implemented; no family-adjusted generalization claim.",
    "Joung OE-vector rebuild from sums302 NOT implemented: 247 measured-axis genes are absent there; no padding or raw-atlas fallback.",
    "Derangements are available as a pure helper but not reported by this driver."
  ],
  "provenance": {
    "manifest_sha256": "7298dfe5258fc1b5bafe5576b7317178fa0d7f9233cfd9041625aa452de4bc48",
    "release_sha256": "81637e9c51d5e364164cda534fd704e7209cd609481dd10d6ef045db1611c78b",
    "product_sha256": "40065fadf89afbca554cc3f88944d9269a05be087e0f92b2e7138738b90ad599",
    "access_sha256": "246c688b7be083fdcd532a716dc3edc90316e9faa0866469e1ca9ea53aae5da8",
    "access_started_at": "2026-09-07T06:01:00Z",
    "source_sha256": "0754d6206fe045744dfd9b21eb37c8e3435152c76b158b9810a84c64a029c0ac",
    "input_sha256": {
      "experiments/303-second-overexpression-library/work/AXIS.txt": "ec1a5719fc5a0bb624e15c13d1a3787a4d89d9c59af31c3dc5e447ef7ecdda20",
      "experiments/303-second-overexpression-library/work/COHORT.tsv": "3d114511339617d16c6c7fde16712c2d5297a52df737f8367966733029a009ea",
      "experiments/303-second-overexpression-library/work/R5-ROWS.tsv": "0c92b72da253e865eb140d0bcad8c41b1b17e315b64dd4c49a457cf46d51ab9b",
      "experiments/303-second-overexpression-library/work/meta303.json": "3515366a01de83eaf3c851f9dede2311f07120599311be0d1bdceb412fa54bd3",
      "experiments/115-gse295214-h1/work/compact-donor/delta.npy": "d0fa6555f9a62ce503d2ad55c8bfda954322657b52dcd48a24e29c330e718c7a",
      "experiments/115-gse295214-h1/work/compact-donor/measured_mask.npy": "da468865f861be0714037c33357b93b9f747351f7478866e6e0dd6dff9d1ebcb",
      "experiments/115-gse295214-h1/work/compact-donor/genes.txt": "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723",
      "experiments/115-gse295214-h1/work/compact-donor/perts.txt": "7060076272d411d4ee00b898cc26eedcb1eaaac5f5e94e2c60f3b45978320747",
      "experiments/279-znf-state-dependence/work/held279.npz": "570b1adf6d55e3b35ae3e3ef25543ec8eef1f2a256836d154154230ec2b51783",
      "experiments/290-remainder-response-features/work/persrc290.npz": "3ce07e88d0aff6330461b350594a4a2202ba8e9af250a93dc1473ef22b19c27c",
      "experiments/279-znf-state-dependence/work/znf279.csv": "f19bfa84a44377671bf66c9f8c868e3a430457065ced3e5ec57f125dc89e58bb",
      "experiments/302-orf-vs-crispri-symmetry/work/symmetry302.csv": "9c15824aed80d984b42a548ab5f8b46c4929adec1752038c6bf1b0bf043ace9f",
      "experiments/303-second-overexpression-library/analyse303.py": "9bf41407a31ff96ec0b99498624f7b20d9aee84dfcc420f132c6fae54b88178a",
      "experiments/303-second-overexpression-library/stats303.py": "9fa088e52523ed22eceaa176eb8d19823db92d6035503864b48fb614ca17cbbe",
      "experiments/303-second-overexpression-library/test_analysis303.py": "ee7a7b4513b8e3375a76bf32c249d2ff9c1ff00cfaf5e61cd08abad755018d97"
    },
    "extraction_code_sha256": {
      "experiments/303-second-overexpression-library/access303.py": "0b1852d4cceb3a5d8b90d3dc6387879baac19b94374f8f0bf28d087e9676ba61",
      "experiments/303-second-overexpression-library/lib303.py": "1ffe619486b1ca9b714eecd10f1fb52d881a2075cc0446223ac23f54d5a9fe84"
    },
    "readout_version": "C7-fixed-axis-uniform-v2",
    "python": "3.12.12",
    "numpy": "2.5.2",
    "scipy": "1.18.0",
    "final_preregistration_sha256": "e92006a21674982d351f1fb17c0ef831577bc08f63da33d831247f66fde3e342"
  }
}
```
