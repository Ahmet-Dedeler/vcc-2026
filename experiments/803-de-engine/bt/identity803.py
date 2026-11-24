"""803 §8 bulk_transform known answer: the identity (rows returned unchanged). Must reproduce E0's runs bit for bit."""


def transform(eff, genes, targets, surface, **args):
    return eff, {"identity": True}
