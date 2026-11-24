"""413: keyword filter shared by the bioRxiv crawlers."""
import re
KW = re.compile(r"perturb-?seq|crop-?seq|crispri|crispra|crispr interference|crispr activation|single-cell crispr|"
                r"pooled crispr|genome-(wide|scale) (crispr|perturb)|perturbation atlas|perturbation screen|"
                r"transcriptom\w* (screen|readout)|l1000|virtual cell", re.I)
