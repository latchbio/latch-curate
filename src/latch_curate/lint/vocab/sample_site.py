sample_site_vocab = {
    "lesional",
    "peri-lesional",
    "normal",
    "blood",
}


def validate_sample_site(val):
    return val in sample_site_vocab
